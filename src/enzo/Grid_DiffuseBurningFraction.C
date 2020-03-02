/***********************************************************************
 /
 /  GRID CLASS (Compute and apply simple burning)
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  based on: ConductHeat and ComputeConductionTimeStep
 /  (by:  David A. Ventimiglia and Brian O'Sheai, Nov, Dec 2009)
 /
 /  PURPOSE:  Computes the new mass density (Rho_56Ni) due to burning
 /  diffusion and constant burning rate and updates the total energy
 /  density.
 /  Updates the gas energy if dual energy folrmalism is on.
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include "DebugMacros.h"
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

#define MEAN_M_12C_16O_MIX_cgs (( M_12C_cgs + M_16O_cgs ) / 2 )
#define MEAN_M_12C_16O_MIX_amu (( M_12C_amu + M_16O_amu ) / 2 )

int rangeLimitInc(float x, float dx, float lowerLimit, float upperLimit,
float *xNew, float *dxUsed);
int lowerLimitInc(float x, float dx, float lowerLimit,
float *xNew, float *dxUsed);
int upperLimitInc(float x, float dx, float upperLimit,
float *xNew, float *dxUsed);

int grid::DiffuseBurnedFraction()
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(NumberOfBaryonFields == 0)
		return SUCCESS;

	switch(BurningDiffusionMethod)
	{
	case -1:
	case 3:
		return SUCCESS;
	case 0:
	case 1:
	case 2:
		break;
	default:
		ENZO_VTHROW(" Unknown diffusion method, %lld, for burned fraction. Valid methods are:\n"
				"   -1    Turn off diffusion. (Different from diffusion rate = reaction rate = 0);\n"
				"    0    7-point 3D stencil;\n"
				"    1    27-point 3D stencil;\n"
				"    2    125-pont 3D stencil.\n",
				BurningDiffusionMethod);
	}

	this->DebugCheck("DiffuseBurnedFraction");

	size_t gridSize = GetGridSize();

	float* F = new float[gridSize];  // Evolving fraction (during the subcycle)
	float* dFdt = new float[gridSize];

	float *Rho, *TE, *GE, *Rho_56Ni;
	MHD_SNIA_GetFields(&Rho, &TE, &GE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &Rho_56Ni, NULL, NULL);

	int Nsub, retcode;
	float dtSubcycle, dtSoFar;
	float f, df, fNew, dfNew;
	float rho, rho_Ni, rho_CO, n_Ni, n_CO, n_all, f_rho, f_mole;
	float X0, X1, Y0, Y1, X1_old, dX1;
	float Q, TE_new, GE_new, Q_for_GE;
	float diffusionRate = BurningDiffusionRateReduced * CellWidth[0][0];
	float reactionRate = BurningReactionRateReduced / CellWidth[0][0];
	float minRhoForBurning = BurningNonDistributedMinDensity;
	double minFractionForDiffusion = BurningMinFractionForDiffusion;

	if(minFractionForDiffusion > 0)
	{
		minFractionForDiffusion *= BurningDiffusionRateReduced / CellWidth[0][0];
		if(minFractionForDiffusion > 0)
			minFractionForDiffusion = powl(minFractionForDiffusion, max(NumberOfBufferZones, NumberOfGhostZones));
	}

	for(size_t i = 0; i < gridSize; i++)
	{
		rho = Rho[i];

		if(rho < tiny_number)
		{
			F[i] = 0.0;
			continue;
		}

		X1 = Rho_56Ni[i] / rho; // Mass fraction of 56Ni
		if(X1 < 0.0)
			X1 = 0.0;
		else if(X1 > 1.0)
			X1 = 1.0;

		X0 = 1 - X1; // Mass fraction of fuel

		Y0 = X0 / MEAN_M_12C_16O_MIX_amu; // Fuel abundance
		Y1 = X1 / M_56Ni_amu; // 56Ni abundance

		F[i] = Y1 / (Y0 + Y1); // Mole fraction of 56Ni
	}

	if(minFractionForDiffusion > 0)
	{
		for(size_t i = 0; i < gridSize; i++)
			if(F[i] < minFractionForDiffusion)
				dFdt[i] = F[i] = 0;
	}

	dtSoFar = 0.0; // cumulative subcycles time, 0 <= dtSoFar <= dtFixed
	Nsub = 0;  // number of subcycles
	// dtSubcycle = timestep of this subcycle
	// dtSoFar = overall timestep taken in heat equation
	// this->dtFixed = the current timestep for the entire level.
	//subcycle begins
	while(dtSoFar < dtFixed)
	{
		if(this->ComputeBurningFractionDiffusionTimeStep(&dtSubcycle) == FAIL)
		{
			ENZO_VFAIL("Grid::DiffuseBurnedFraction: Error in ComputeBurningTimeStep. %e %e\n", dtFixed, dtSubcycle);
		}
		// Make sure we don't extend past dtFixed
		dtSubcycle = min(dtSubcycle, dtFixed - dtSoFar);

		//Compute the diffusion term
		if(this->ComputeLaplacian(dFdt, F, BurningDiffusionMethod) == FAIL)
			ENZO_VFAIL("Grid::DiffuseBurnedFraction: Error in ComputeLaplacian. dtSub,dtFixed: %e %e\n", dtFixed,
						dtSubcycle);

		arr_ax(dFdt, gridSize, diffusionRate);

		if(reactionRate)
		{
			// Add the source term
			for(size_t i = 0; i < gridSize; i++)
				if(BurningReactionBurnedFractionLimitLo < F[i] && F[i] < BurningReactionBurnedFractionLimitHi)
					dFdt[i] += reactionRate;
		}

		// Evolve F
		for(size_t i = 0; i < gridSize; i++)
		{
			int retcode = rangeLimitInc(F[i], dFdt[i] * dtSubcycle, 0, 1, &f, &df);
			F[i] = f;
		}

		if(minFractionForDiffusion > 0)
		{
			for(size_t i = 0; i < gridSize; i++)
				if(F[i] < minFractionForDiffusion)
					dFdt[i] = F[i] = 0;
		}

		dtSoFar += dtSubcycle;
		Nsub++;
	} // while(dtSoFar < dtFixed)
	  // end of subcycle

	// Update the burned density and internal energy using the new burned fraction, F.
	if(debug1)
		printf("Grid::DiffuseBurnedFraction: Nsubcycles = %"ISYM"\n", Nsub);

	float totQp=0, totQm=0;
	int skip1 = 0, skip2 = 0, skip3 = 0, skip4 = 0;
	for(size_t i = 0; i < gridSize; i++)
	{
		rho = Rho[i];
		if(rho < tiny_number)
		{
			Rho_56Ni[i] = 0;
			continue;
		}

		X1_old = Rho_56Ni[i] / rho;

		f_mole = F[i];
		X1 = f_mole * M_56Ni_amu; // intermediate result
		X1 /= ( MEAN_M_12C_16O_MIX_amu * (1 - f_mole) + X1); // New mass fraction of 56Ni
		// Calculate the new mass fraction of 56Ni making sure it stays between 0 and 1.
		retcode = rangeLimitInc(X1_old, X1 - X1_old, 0, 1, &X1, &dX1);

		Rho_56Ni[i] = X1 * rho; // Update the 56Ni mass density.

		Q = dX1 * BurningEnergyRelease; // Enrgy release per gram product


		if(rho < minRhoForBurning)
		{
			skip1++;
			continue;
		}

		Q = dX1 * BurningEnergyRelease; // Enrgy release per gram product

		if(Q == 0)
		{
			skip2++;
			continue;
		}
		if(Q < 0 && !AllowUnburning)
		{
			skip3++;
			continue;
		}
		skip4++;
		if(Q>0) totQp += Q;
		if(Q<0) totQm += Q;
//		if(Q > 0)
//			TRACEGF("  TE=%e  Q=%e  dX1=%e  Q0=%e", TE[i], Q, dX1, BurningEnergyRelease);
		TE[i] += Q;
		if(GE)
			GE[i] += Q;
	}
	TRACEGF("  totQp=%e  totQm=%e", totQp, totQm);
	TRACEGF(" skippers:  %lld  (%e) %lld  %lld  (%lld)  %lld  (%e)", skip1, minRhoForBurning, skip2, skip3, AllowUnburning, skip4, BurningEnergyRelease);
	{
		int i;
		if(GE)
		{
			if(0 <= (i = findmaxlte(GE, gridSize, 0)))
			{
				TRACEGF("Non-positive internal energy, %e, at index %lld", GE[i]);
			}
		}

		if(0 <= (i = findmaxlte(TE, gridSize, 0)))
		{
			TRACEGF("Non-positive total energy, %e, at index %lld", TE[i]);
		}
	}

	delete[] dFdt;
	delete[] F;

	if(InitialBurnedRegionSustain)
		MHDSustainInitialBurnedRegionGrid();

	return SUCCESS;

}

//#endif
/***********************************************************************
 /
 /  GRID CLASS (Compute burning diffusion time step)
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  based on: ConductHeat and ComputeConductionTimeStep
 /  (by:  David A. Ventimiglia and Brian O'Sheai, Nov, Dec 2009)
 /
 /  PURPOSE:  Calculates a safe time step for the subcycle evolving
 /    the burned fraction diffusion equation.
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

// Function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
float *TemperatureUnits, float *TimeUnits,
float *VelocityUnits, double *MassUnits, FLOAT Time);

// Member functions
float grid::ComputeBurningFractionDiffusionTimeStep(float* dt)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	this->DebugCheck("Grid::ComputeBurningFractionDiffusionTimeStep");

	// Some locals
	float dt1;
	float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
	float VelocityUnits = 1.0, TimeUnits = 1.0, aUnits = 1.0;
	double MassUnits = 1.0;
	float light_cross_time;
	double all_units;

	FLOAT dx = CellWidth[0][0];					// For some reason dx==0 here. [BH]
	for(int dim = 1; dim < GridRank; dim++)
		dx = min(dx, CellWidth[dim][0]);

	if(GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL)
		ENZO_FAIL("Error in GetUnits.");

	//TODO verify use of units in the following few statements
	//  dt1 = BurningDiffusionCourantSafetyFactor * POW( dx * LengthUnits , 2.0)
	//       / ( diffusionRate / TimeUnits );
	dt1 = BurningDiffusionCourantSafetyFactor * dx * LengthUnits / (BurningDiffusionRateReduced / TimeUnits);

//printf("BurningDiffusionCourantSafetyFactor=%e\n", BurningDiffusionCourantSafetyFactor);
//printf("dx=%e\n", dx);
//printf("LengthUnits=%e\n", LengthUnits);
//printf("BurningDiffusionRate=%e\n", BurningDiffusionRate);
//printf("TimeUnits=%e\n", TimeUnits);
//printf("dt=%e, t=%e\n", dt, dtFixed );
//ENZO_FAIL("debug");

//  if (SpeedOfLightTimeStepLimit) {
//    light_cross_time = dx * VelocityUnits / clight;
//    dt1 = max(dt1, light_cross_time);
//  }

	*dt = dt1;

	return SUCCESS;
}

/***********************************************************************
 /
 /  Range limit increment
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx stays in
 /	the specified range [lowerLimit, upperLimit].
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxUsed=dx and xNew=x+dx
 /   -1 - 'dx' was too negative causing 'x' to undeshoot
 /	Output values: dxUsed=lowerLimit-x and xNew=lowerLimit
 /    1 - 'dx' was too positiv) causing 'x' to overshot
 /	Output values: dxUsed=upperLimit-x and xNew=upperLimit
 /
 ************************************************************************/
int rangeLimitInc(float x, float dx, float lowerLimit, float upperLimit,
float *xNew, float *dxUsed)
{
	float y = x + dx;

	if(y < lowerLimit)
	{
		*dxUsed = (*xNew = lowerLimit) - x;
		return -1;
	}

	if(y > upperLimit)
	{
		*dxUsed = (*xNew = upperLimit) - x;
		return 1;
	}

	*dxUsed = dx;
	*xNew = y;
	return 0;
}

/***********************************************************************
 /
 /  Lower limit increment
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx doesn't go
 /	below lowerLimit.
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxUsed=dx and xNew=x+dx
 /   -1 - 'dx' was too negative causing 'x' to undeshoot
 /	Output values: dxUsed=lowerLimit-x and xNew=lowerLimit
 /
 ************************************************************************/
int lowerLimitInc(float x, float dx, float lowerLimit,
float *xNew, float *dxUsed)
{
	float y = x + dx;

	if(y < lowerLimit)
	{
		*dxUsed = (*xNew = lowerLimit) - x;
		return -1;
	}

	*dxUsed = dx;
	*xNew = y;
	return 0;
}

/***********************************************************************
 /
 /  Upper limit increment
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx doesn't go
 /	above upperLimit.
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxUsed=dx and xNew=x+dx
 /    1 - 'dx' was too positive causing 'x' to overshoot
 /	Output values: dxUsed=upperLimit-x and xNew=upperLimit
 /
 ************************************************************************/
int upperLimitInc(float x, float dx, float upperLimit,
float *xNew, float *dxUsed)
{
	float y = x + dx;

	if(y > upperLimit)
	{
		*dxUsed = (*xNew = upperLimit) - x;
		return 1;
	}

	*dxUsed = dx;
	*xNew = y;
	return 0;
}

