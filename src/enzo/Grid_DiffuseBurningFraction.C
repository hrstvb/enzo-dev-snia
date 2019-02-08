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

int rangeCutLimiter(float x, float dx, float lowCutoff, float highCutoff,
float *xNew, float *dxNew);
int lowCutLimiter(float x, float dx, float lowCutoff,
float *xNew, float *dxNew);
int highCutLimiter(float x, float dx, float highCutoff,
float *xNew, float *dxNew);

int grid::DiffuseBurnedFraction()
{

	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(NumberOfBaryonFields == 0)
		return SUCCESS;

	this->DebugCheck("DiffuseBurnedFraction");

	int gridSize = GetGridSize();

	float* F = new float[gridSize];  // Evolving fraction (during the subcycle)
	float* dFdt = new float[gridSize];

	int useGE = DualEnergyFormalism && (HydroMethod != Zeus_Hydro);

	int RhoNum = FindField(Density, FieldType, NumberOfBaryonFields);
	int TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
	int GENum = (useGE) ? FindField(InternalEnergy, FieldType, NumberOfBaryonFields) : 0;
	int Rho_56NiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);

	if(RhoNum < 0)
		ENZO_FAIL("Grid::DiffuseBurnedFraction: Cannot find density.");
	if(TENum < 0)
		ENZO_FAIL("Grid::DiffuseBurnedFraction: Cannot find total energy.");
	if(useGE && GENum < 0)
		ENZO_FAIL("Grid::DiffuseBurnedFraction: Cannot find gas energy.");
	if(Rho_56NiNum < 0)
		ENZO_FAIL("Grid::DiffuseBurnedFraction: Cannot find Density_56Ni.");

	float *Rho = BaryonField[RhoNum];
	float *TE = BaryonField[TENum];
	float *GE = (useGE) ? BaryonField[GENum] : NULL;
	float *Rho_56Ni = BaryonField[Rho_56NiNum];

	int idx, i, j, k, Nsub, retcode;
	int fractionMode;
	float dtSubcycle, dtSoFar;
	float f, df, fNew, dfNew;
	float rho, rho_Ni, rho_CO, n_Ni, n_CO, n_all, f_rho, f_mole;
	float X0, X1, Y0, Y1, X1_old, dX1;
	float Q, TE_new, Q_new_TE, GE_new, Q_new_GE;
	float diffusionRate = BurningDiffusionRateReduced * CellWidth[0][0];
	float reactionRate = BurningReactionRateReduced / CellWidth[0][0];

	// dtSubcycle = timestep of this subcycle
	// dtSoFar = overall timestep taken in heat equation
	// this->dtFixed = the current timestep for the entire level.

	for(i = 0; i < gridSize; i++)
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

	dtSoFar = 0.0;

	Nsub = 0;  // number of subcycles

	//subcycle begins
	while(dtSoFar < dtFixed)
	{

		if(this->ComputeBurningFractionDiffusionTimeStep(&dtSubcycle) == FAIL)
		{
			ENZO_VFAIL("Grid::DiffuseBurnedFraction: Error in ComputeBurningTimeStep. %e %e\n", dtFixed, dtSubcycle);
		}

		// make sure we don't extend past dtFixed
		dtSubcycle = min(dtSubcycle, dtFixed - dtSoFar);

		if(this->ComputeLaplacian(dFdt, F, BurningDiffusionMethod) == FAIL)
			ENZO_VFAIL("Grid::DiffuseBurnedFraction: Error in ComputeLaplacian. dtSub,dtFixed: %e %e\n", dtFixed,
						dtSubcycle);

		arr_ax(dFdt, gridSize, diffusionRate);

		if(reactionRate)
		{
			for(i = 0; i < gridSize; i++)
			{
//			dFdt[i] *= diffusionRate;

				// Add reaction
//			dFdt[i] +=
//					(BurningReactionBurnedFractionLimitLo < F[i] && F[i] < BurningReactionBurnedFractionLimitHi) ?
//							reactionRate : 0;
				if(BurningReactionBurnedFractionLimitLo < F[i] && F[i] < BurningReactionBurnedFractionLimitHi)
					dFdt[i] += reactionRate;

				// Evolve F
				int retcode = rangeCutLimiter(F[i], dFdt[i] * dtSubcycle, 0, 1, &f, &df);
//TODO: change error to warning
//      if ( retcode ) {
//	    ENZO_VFAIL("Grid::DiffuseBurnedFraction: bad burned fraction=%g dFdt=%g i=%d  dtFixed,dtSub: %e, %e\n",
//		       F[i], dFdt[i], i, dtFixed, dtSubcycle)
//      }
				F[i] = f;

			} // end for i
		}
		dtSoFar += dtSubcycle;
		Nsub++;

	} // while(dtSoFar < dtFixed)
// end of subcycle

	if(debug1)
		printf("Grid::DiffuseBurnedFraction: Nsubcycles = %"ISYM"\n", Nsub);

//*LevelArray[0]->GridData->*/ExtraFunction("[BH]:>>>>>>>> SolveMHDLi BEFORE LAPL", 0, NULL, NULL); //[BH]EF
	for(i = 0; i < gridSize; i++)
	{
		rho = Rho[i];

		if(rho < tiny_number)
		{
			//This is an empty cell (rho=0).
			//TODO: Decide betweeni: (a) Set TE=GE=0; and (b) Do nothing, i.e. leave TE and GE unchanged.
			continue;
		}

		X1_old = Rho_56Ni[i] / rho;

		f_mole = F[i];
		X1 = f_mole * M_56Ni_amu; // intermediate result
		X1 /= ( MEAN_M_12C_16O_MIX_amu * (1 - f_mole) + X1); // New mass fraction of 56Ni
		// Calculate the new mass fraction of 56Ni making sure it stays between 0 and 1.
		retcode = rangeCutLimiter(X1_old, X1 - X1_old, 0, 1, &X1, &dX1);

		Rho_56Ni[i] = X1 * rho; //Update the 56Ni mass density field with the new value.
		Q = dX1 * BurningEnergyRelease; //enrgy release per gram product

		if(Q > 0 || Q < 0 && AllowUnburning) // added 11/19/2015 [BH]
		{

			// In case Q is negative make sure energy doesn't get negative:
			retcode = lowCutLimiter(TE[i], Q, 0, &TE_new, &Q_new_TE);
			//TODO: change error to warning
			//      if ( retcode < 0 ) {
			//        ENZO_VFAIL("Grid::DiffuseBurnedFraction: Bad total energy=%g dFdt=%g i=%d  dtFixed,dtSub: %e, %e\n",
			//             dE[i], df, i, dtFixed, dtSubcycle)
			//      }

			if(useGE)
			{
				retcode = lowCutLimiter(GE[i], Q, 0, &GE_new, &Q_new_GE);
				//TODO: change error to warning
				//      if ( retcode < 0 ) {
				//          ENZO_VFAIL("Grid::DiffuseBurnedFraction: Bad gas energy=%g dFdt=%g i=%d  dtFixed,dtSub: %e, %e\n",
				//                    GE[i], df, i, dtFixed, dtSubcycle)
				//      }

				//Q = max( Q_new_TE, Q_new_GE ); replaced 11/19/2015 [BH]
				Q_new_TE = max(Q_new_TE, Q_new_GE); // max(...) because we worry only about negative values
				//GE[i] += Q; replaced 11/19/2015 [BH]
				GE[i] += Q_new_TE;
			}

			//TE[i] += Q_new_TE; replaced 11/19/2015 [BH]
			TE[i] += Q_new_TE;
		}
		else
		{
		}
	}
//*LevelArray[0]->GridData->*/ExtraFunction("[BH]:>>>>>>>> SolveMHDLi AFTER LAPL", 0, NULL, NULL); //[BH]EF

	delete[] dFdt;
	delete[] F;

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

	FLOAT dx = CellWidth[0][0]; // For some reason dx==0 here. [BH]
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
 /  Range-cut limiter
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx stays in
 /	the specified range [lowCutoff, highCutoff].
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxNew=dx and xNew=x+dx
 /   -1 - 'dx' was too big (too negative) causing 'x' to undeshoot
 /	Output values: dxNew=lowCutoff-x and xNew=lowCutoff
 /    1 - 'dx' was too big (too positive) causing 'x' to overshot
 /	Output values: dxNew=highCutoff-x and xNew=highCutoff
 /
 ************************************************************************/
int rangeCutLimiter(float x, float dx, float lowCutoff, float highCutoff,
float *xNew, float *dxNew)
{
	float y = x + dx;

	if(y < lowCutoff)
	{
		*dxNew = (*xNew = lowCutoff) - x;
		return -1;
	}

	if(y > highCutoff)
	{
		*dxNew = (*xNew = highCutoff) - x;
		return 1;
	}

	*dxNew = dx;
	*xNew = y;
	return 0;
}

/***********************************************************************
 /
 /  Low-cut limiter
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx doesn't go
 /	below lowCutoff.
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxNew=dx and xNew=x+dx
 /   -1 - 'dx' was too big (too negative) causing 'x' to undeshoot
 /	Output values: dxNew=lowCutoff-x and xNew=lowCutoff
 /
 ************************************************************************/
int lowCutLimiter(float x, float dx, float lowCutoff,
float *xNew, float *dxNew)
{
	float y = x + dx;

	if(y < lowCutoff)
	{
		*dxNew = (*xNew = lowCutoff) - x;
		return -1;
	}

	*dxNew = dx;
	*xNew = y;
	return 0;
}

/***********************************************************************
 /
 /  High-cut limiter
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Limits the growth 'dx' so that xNew=x+dx doesn't go
 /	above highCutoff.
 /
 /  RETURNS:
 /    0 - no adjustment of 'dx'
 /	Output values: dxNew=dx and xNew=x+dx
 /    1 - 'dx' was too big (too positive) causing 'x' to overshot
 /	Output values: dxNew=highCutoff-x and xNew=highCutoff
 /
 ************************************************************************/
int highCutLimiter(float x, float dx, float highCutoff,
float *xNew, float *dxNew)
{
	float y = x + dx;

	if(y < highCutoff)
	{
		*dxNew = (*xNew = highCutoff) - x;
		return 1;
	}

	*dxNew = dx;
	*xNew = y;
	return 0;
}

