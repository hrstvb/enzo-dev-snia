/***********************************************************************
 /
 /  GRID CLASS (INITIALIZE THE MHD BLAST GRID)
 /
 /  written by: David Collins
 /  date:       2004-2013
 /
 /  modified1:  Boyan Hristov
 /  date:       2015
 /  		Added perturbation method 33 creating a plane interface
 /		between zones A and B perturbed by a sine offset.
 /
 /  PURPOSE:  See MHDBlastInitialize.C for parameters
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <cctype> //isspace

#include "myenzoutils.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

#include "DebugMacros.h"
#include "MHDInitialProfile.h"

//#ifndef SRC_ENZO_INITIALPROFILE_H_
//#define SRC_ENZO_INITIALPROFILE_H_
//
//#define LINE_MAX_LENGTH (512)
//#define CMP_A_B(a, b) ((b>a)-(b<a))
//#define SIGN_A(a) (CMP_A_B(0,a))
//
//#define PROFILE_ASC_SORT (1)
//#define PROFILE_DESC_SORT (-1)
//#define PROFILE_UNSORTED (0)
//
//#define PROFILE_EMPTY_LINE (0)
//#define PROFILE_COMMENT_LINE (1)
//#define PROFILE_COL_NAMES_LINE (2)
//#define PROFILE_DATA_LINE (3)
//#define PROFILE_TIME_LINE (4)

//using namespace std;

int MakeFieldConservative(int field);
int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z);
float SphericalGravityGetAt(FLOAT r);
size_t SphericalGravityComputeBinIndex(FLOAT r);
void WriteInitialProfile(char* name, FLOAT* RR, FLOAT* RHO, FLOAT* GG, FLOAT* PP, FLOAT* UU, size_t n, FLOAT K,
FLOAT gamma);

int grid::MHD_SNIA_GetFields(float** densityField, float** totalEnergyField, float** internalEnergyField,
float** vxField, float** vyField, float** vzField, float** vFields,
float** BxField, float** ByField, float** BzField, float** BFields,
float** burnedDensityField,
float** phiField,
float** gravPotentialField)
{
	bool hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	int totENum, rhoNum, vxNum, vyNum, vzNum;
	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	int phiNum = -1, gravPotentialNum = -1;

	if(IdentifyPhysicalQuantities(rhoNum, gasENum, vxNum, vyNum, vzNum, totENum, BxNum, ByNum, BzNum, phiNum) == FAIL)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in IdentifyPhysicalQuantities for UseMHD==true.")

	if(densityField)
		*densityField = BaryonField[rhoNum];

	if(totalEnergyField)
		*totalEnergyField = BaryonField[totENum];

	if(internalEnergyField)
		*internalEnergyField = (hasGasEField) ? BaryonField[gasENum] : NULL;

	if(vxField)
		*vxField = BaryonField[vxNum];
	if(vyField)
		*vyField = (MaxVelocityIndex > 1) ? BaryonField[vyNum] : NULL;
	if(vzField)
		*vzField = (MaxVelocityIndex > 2) ? BaryonField[vzNum] : NULL;
	if(vFields)
	{
		vFields[0] = BaryonField[vxNum];
		vFields[1] = (MaxVelocityIndex > 1) ? BaryonField[vyNum] : NULL;
		vFields[2] = (MaxVelocityIndex > 2) ? BaryonField[vzNum] : NULL;
	}

	if(UseMHD)
	{
		if(BxField)
			*BxField = BaryonField[BxNum];
		if(ByField)
			*ByField = BaryonField[ByNum];
		if(BzField)
			*BzField = BaryonField[BzNum];

		if(BFields)
		{
			BFields[0] = BaryonField[BxNum];
			BFields[1] = BaryonField[ByNum];
			BFields[2] = BaryonField[BzNum];
		}
	}

	if(burnedDensityField && UseBurning)
	{
		rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
		if(rhoNiNum < 0)
			ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
		*burnedDensityField = BaryonField[rhoNiNum];
	}

	if(phiField)
		*phiField = (phiNum >= 0) ? BaryonField[phiNum] : NULL;

	if(gravPotentialField && WritePotential)
	{
		gravPotentialNum = FindField(GravPotential, FieldType, NumberOfBaryonFields);
		if(gravPotentialNum < 0)
			ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for GravPotential.")
		*gravPotentialField = BaryonField[gravPotentialNum];
	}

	return SUCCESS;
}

int grid::MHDProfileInitializeGrid(MHDInitialProfile* p,
float burningTemperature,
float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool usingVectorPotential)
{
//	TRACEGF("INITIALIZING GRID PHASE 1 ...")
	if(GridRank != 3)
		ENZO_FAIL("MHDProfileInitializeGrid is implemented for 3D only.")

	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	if(CellWidth[0][0] <= 0)
		PrepareGridDerivedQuantities();

// Assign fieldType numbers using constants from typedefs.h,
// as well as count the number of fields.
	NumberOfBaryonFields = 0;
	FieldType[NumberOfBaryonFields++] = Density;
	if(EquationOfState == 0)
		FieldType[NumberOfBaryonFields++] = TotalEnergy;
	if(hasGasEField)
		FieldType[NumberOfBaryonFields++] = InternalEnergy;
	FieldType[NumberOfBaryonFields++] = Velocity1;
	FieldType[NumberOfBaryonFields++] = Velocity2;
	FieldType[NumberOfBaryonFields++] = Velocity3;
	if(WritePotential)
		FieldType[NumberOfBaryonFields++] = GravPotential;
	if(UseBurning)
		FieldType[NumberOfBaryonFields++] = Density_56Ni;   //[BH]
	if(UseMHD)
	{
		FieldType[NumberOfBaryonFields++] = Bfield1;
		FieldType[NumberOfBaryonFields++] = Bfield2;
		FieldType[NumberOfBaryonFields++] = Bfield3;
		if(HydroMethod == MHD_RK)
			FieldType[NumberOfBaryonFields++] = PhiField;
	}

	if(HydroMethod == MHD_RK || usingVectorPotential)
	{
		// Allow the ElectricField and MagneticField to
		// be created temporarily for the purpose of
		// initializing with vector potential.
		// Free these fields in MHDProfileInitializeGrid2.
		UseMHDCT = TRUE;
		MHD_SetupDims();
	}

	if(ProcessorNumber != MyProcessorNumber || p == NULL)
	{
		//p==NULL is used with ParallelGridIO.
//		TRACEGF("... DONE INITIALIZING GRID W/O FIELDS PHASE 1.")
		return SUCCESS;
	}

	TRACEGF("INITIALIZING GRID PHASE 1, FIELDS.")

	this->AllocateGrids();

//	int totENum, rhoNum, vxNum, vyNum, vzNum;
//	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField, NULL, &rhoNiField, NULL, NULL);

	int debugnanflag = 0; //[BH]
	float gammaMinusOne = Gamma - 1;

	for(int k = 0; k < GridDimension[2]; k++)
	{
		FLOAT z = CELLCENTER(2, k);
		FLOAT zz = square(z - SphericalGravityCenter[2]);
		for(int j = 0; j < GridDimension[1]; j++)
		{
			FLOAT y = CELLCENTER(1, j);
			FLOAT yy_zz = square(y - SphericalGravityCenter[1]) + zz;
			for(int i = 0; i <= GridDimension[0]; i++)
			{
				int index = ELT(i, j, k);
				float vx = vxField[index] = 0 + GridVelocity[0];
				float vy = vyField[index] = 0 + GridVelocity[1];
				float vz = vzField[index] = 0 + GridVelocity[2];
				float Bx, By, Bz;

				if(1) //!strcmp(p->profileType, "RADIAL"))
				{
					FLOAT x = CELLCENTER(0, i);
					FLOAT r = sqrt(square(x - SphericalGravityCenter[0]) + yy_zz);

					float rho, T = 0;
					int retcode = p->interpolateDensity(&rho, r); //g/cm**3

					//retcode = profileInterpolate(&T, temperatureData, r, radiusData, p.nRows, radiusSO); //K
					bool isBurned = (r < burnedRadius);				// || (T > 0 && T > burningTemperature));

					float gasE = 0;
					if(p->internalEnergyData)
					{
						retcode = p->interpolateInternalEnergy(&gasE, r); //g/cm**3
					}
					else
					{
						gasE = EOSPolytropicFactor * POW(rho, gammaMinusOne) / gammaMinusOne;
					}

//					if(UseBurning && isBurned && burnedRadiusPE)
//					{
//						rho *= gasE / (gasE + BurningEnergyRelease);
//						gasE += BurningEnergyRelease;
//					}

					float totE = gasE;
					totE += 0.5 * (square(vx) * square(vy) + square(vz));
					if(BxField)
					{
						MHDProfileInitExactB(&Bx, &By, &Bz, x, y, z);
						totE += 0.5 * (square(Bx) + square(By) + square(Bz)) / rho;
						BxField[index] = Bx;
						ByField[index] = By;
						BzField[index] = Bz;
					}

					rhoField[index] = rho;
					if(hasGasEField)
						gasEField[index] = gasE;
					totEField[index] = totE;
					vxField[index] = vx;
					vyField[index] = vy;
					vzField[index] = vz;

					if(UseBurning)
						rhoNiField[index] = (isBurned) ? rho : 0;

					if(debug + 1 && MyProcessorNumber == ROOT_PROCESSOR)
						if(j == (GridDimension[1] / 2) && k == (GridDimension[2] / 2) && i >= GridDimension[0] / 2)
							TRACEF("i,j,k=%03d,%04d,%04d, x,y,z=(%4f,%4f,%4f), r=%4f, rho=%e, U=%e, T=%1f, burned=%d, K=%e, gamma-1=%f",
									i, j, k, x * 1e-5, y * 1e-5, z * 1e-5, r * 1e-5, rho, gasE, T * 1e-9, isBurned,
									EOSPolytropicFactor, gammaMinusOne);
				}
			} //end baryonfield initialize
		}
	}

// Boiler plate code:
//  if(DualEnergyFormalism )
//    for(index=0;index<size;index++)
//    BaryonField[ gesENum ][index] =
//       BaryonField[ totENum ][index]
//      - 0.5*(BaryonField[ vNum[0] ][index]*BaryonField[ vNum[0] ][index] +
//	     BaryonField[ vNum[1] ][index]*BaryonField[ vNum[1] ][index] +
//	     BaryonField[ vNum[2] ][index]*BaryonField[ vNum[2] ][index])
//      - 0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index] +
//             BaryonField[ByNum][index]*BaryonField[ByNum][index] +
//             BaryonField[BzNum][index]*BaryonField[BzNum][index])/BaryonField[ rhoNum ][index];

	if(debugnanflag)
		printf("nans intialized.\n"); //[BH]
//	else
//		printf("Initialized grid, phase 1.\n");
	TRACEGF("... DONE INITIALIZING GRID FULL PHASE 1.")

	return SUCCESS;
}

/*
 * Keeps the burned fraction equal to 1 for the initial burned raadius.
 */
int grid::MHDMaintaindInitialBurnedRegionGrid()
{
	if(GridRank != 3)
		ENZO_FAIL("MHDProfileInitializeGridBurnedFraction is implemented for 3D only.")

	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(!UseBurning)
		return SUCCESS;

	int rhoNum = FindField(Density, FieldType, NumberOfBaryonFields);
	if(rhoNum < 0)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
	int rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
	if(rhoNiNum < 0)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")

	float* rhoField = BaryonField[rhoNum];
	float* rhoNiField = BaryonField[rhoNiNum];

	FLOAT xyz[MAX_DIMENSION];
	size_t lindex[MAX_DIMENSION];
	arr_set(xyz, MAX_DIMENSION, 0);
	arr_set(xyz, GridRank, InitialBurnedRadius);
	arr_xpy(xyz, SphericalGravityCenter, GridRank);
	get_ijk_index(lindex, xyz);
	for(int dim = 0; dim < GridRank; dim++)
	{
		if(lindex[dim] < 0)
			lindex[dim] = 0;
		else if(lindex[dim] >= GridDimension[dim])
			lindex[dim] >= GridDimension[dim] - 1;
	}

	size_t rindex[MAX_DIMENSION];
	arr_set(xyz, MAX_DIMENSION, 0);
	arr_set(xyz, GridRank, -InitialBurnedRadius);
	arr_xpy(xyz, SphericalGravityCenter, GridRank);
	get_ijk_index(rindex, xyz);
	arr_xpa(rindex, MAX_DIMENSION, 1);
	for(int dim = 0; dim < GridRank; dim++)
	{
		if(rindex[dim] < 0)
			rindex[dim] = 0;
		else if(rindex[dim] >= GridDimension[dim])
			rindex[dim] >= GridDimension[dim] - 1;
	}

	switch(GridRank)
	{
	case 1:
		for(int i = lindex[0]; i <= rindex[0]; i++)
		{
			FLOAT r = sqrt(square(CELLCENTER(0, i) - SphericalGravityCenter[0]));

			if(r < InitialBurnedRadius)
			{
				int index = i;
				rhoNiField[index] = rhoField[index];
			}
		}
		break;
	case 2:
		for(int j = lindex[1]; j < rindex[1]; j++)
		{
			FLOAT yy = square(CELLCENTER(1, j) - SphericalGravityCenter[1]);
			for(int i = lindex[0]; i <= rindex[0]; i++)
			{
				FLOAT r = sqrt(square(CELLCENTER(0, i) - SphericalGravityCenter[0]) + yy);

				if(r < InitialBurnedRadius)
				{
					int index = ELT(i, j, 0);
					rhoNiField[index] = rhoField[index];
				}
			}
		}
		break;
	case 3:
		for(int k = lindex[2]; k <= rindex[2]; k++)
		{
			FLOAT zz = square(CELLCENTER(2, k) - SphericalGravityCenter[2]);
			for(int j = lindex[1]; j < rindex[1]; j++)
			{
				FLOAT yy_zz = square(CELLCENTER(1, j) - SphericalGravityCenter[1]) + zz;
				for(int i = lindex[0]; i <= rindex[0]; i++)
				{
					FLOAT r = sqrt(square(CELLCENTER(0, i) - SphericalGravityCenter[0]) + yy_zz);

					if(r < InitialBurnedRadius)
					{
						int index = ELT(i, j, k);
						rhoNiField[index] = rhoField[index];
					}
				}
			}
		}
		break;
	default:
		ENZO_VFAIL("Bad grid rank: %lld", GridRank)
		;
	}

	return SUCCESS;
}

/*
 * It resets the total energy with the updated magnetic field.
 * As a side effect it frees the E&M fields for hydro method MHD_RK.
 * This function is to be used if and  after magnetic vector
 * potential has been projected to parents and the curl taken.
 * The first initializer still calculates the total energy,
 * as well as possible, because it might influence the
 * start up refinement.
 */
int grid::MHDProfileInitializeGrid2(MHDInitialProfile* p,
float burningTemperature,
float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool usingVectorPotential)
{

	if(ProcessorNumber != MyProcessorNumber || p == NULL)
		return SUCCESS;

	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

//	int radiusSO = p->colSortingOrders[profileFindColIndex(radiusColumnName, p)];
//	float* radiusData = profileFindCol(radiusColumnName, p);
//	float* densityData = profileFindCol(densityColumnName, p);
//	float* gasEData = profileFindCol(InternalEnergyColumnName, p);

	double rho, rho_c;
	int retcode; // = profileInterpolate(&rho, densityData, r, radiusData, p->nRows, 1); //g/cm**3

//	int totENum, rhoNum, vxNum, vyNum, vzNum;
//	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	float *totEField, *rhoField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;
	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField, NULL, &rhoNiField, NULL, NULL);

	float gammaMinusOne = Gamma - 1;
	size_t gridSize = GetGridSize();
	const float* totEField_end = totEField + gridSize;

	float* RHO = rhoField;
	float* GE = gasEField;
	float* VX = vxField;
	float* VY = vyField;
	float* VZ = vzField;
	float* BX = BxField;
	float* BY = ByField;
	float* BZ = BzField;

	if(InitRadialPressureFromCentral)
	{
		// Calculate the pressure for hydrostatic equilibrium with the exisiting gravity.
//		float P_c = 0;
//		float rho_c = 0;
//		FLOAT xyz[MAX_DIMENSION];
//		arr_cpy(xyz, SphericalGravityCenter, MAX_DIMENSION);
//
//		if(InitRadialPressureFromCentral > 0)
//		{
//			P_c = InitRadialPressureFromCentral;
//			if(MyProcessorNumber == ROOT_PROCESSOR)
//				printf("P_c = %e (parameter)", P_c);
//		}
//		else
//		{
////			size_t i_c = getCellIndex(xyz);
////			float rho_c = rhoField[i_c];
//			FLOAT rho_c;
//			retcode = p->interpolateDensity(&rho_c, 0); //g/cm**3
//			P_c = EOSPolytropicFactor * POW(rho_c, Gamma);
//			if(MyProcessorNumber == ROOT_PROCESSOR)
//				printf("P_c = %e (polytropic, = K*rho**gamma, K=%e, rho=%e, gamma=%5.3f)", P_c, EOSPolytropicFactor,
//						rho_c, gamma);
//		}
//
		FLOAT* P_r = new float[SphericalGravityActualNumberOfBins];
		FLOAT* U_r = new float[SphericalGravityActualNumberOfBins];
		FLOAT* rr = NULL;
		FLOAT* dd = NULL;
		FLOAT* gg = NULL;
//		if(ProcessorNumber == ROOT_PROCESSOR && ID == 0 && CellWidth[0][0] == TopGridDx[0])
//		{
//			rr = new FLOAT[SphericalGravityActualNumberOfBins];
//			gg = new FLOAT[SphericalGravityActualNumberOfBins];
//			dd = new FLOAT[SphericalGravityActualNumberOfBins];
//		}
//
		FLOAT P2, P1, P, r2, r1, r, U2, U1, U; // r2 < r1 < r
		FLOAT b, a, bb, aa, g, rho, g_rho; // b:=r-r2 > a:=r-r1
//
//		r2 = 0; // Start from the center.
//		P_r[0] = P2 = 0; // Offfset the pressure later.
//
//		//xyz[0] = r2 + SphericalGravityCenter[0];
//		//size_t index = getCellIndex(xyz);
//		retcode = p->interpolateDensity(&rho, r2); //g/cm**3
//		g_rho = g = SphericalGravityGetAt(r2);
////		TRACEF("r[0], P[0], g[0]rho[0], g[0], rho[0] = %e, %e, %e = %e * %e", r2, P2, g_rho, g, rho);
//		if(rr)
//		{
//			rr[0] = r2;
//			gg[0] = g;
//			dd[0] = rho;
//		}
//
//		r1 = SphericalGravityBinLeftEdges[1];
////		xyz[0] = r1 + SphericalGravityCenter[0];
////		index = getCellIndex(xyz);
////		rho1 = rhoField[index];
//		retcode = p->interpolateDensity(&rho, r1); //g/cm**3
//		g = SphericalGravityGetAt(r1);
//		g_rho = g * rho;
//		if(rr)
//		{
//			rr[1] = r1;
//			gg[1] = g;
//			dd[1] = rho;
//		}
//
//		a = r1 - r;
//		P1 = P2 - 0.5 * a * g * g_rho;
//		P1 -= 2. / 3. * M_PI * SphericalGravityConstant * a * a * rho * rho;
//
//		P_r[1] = P1;
////		TRACEF("r[1], P[1], dP[1], g, rho = %e, %e, %e = %e * %e", r1, P1, P1 - P2, g, rho);
//
//		for(size_t i = 2; i < SphericalGravityActualNumberOfBins; i++)
//		{
//			r = SphericalGravityBinLeftEdges[i];
//			b = r - r2;
//			a = r - r1;
//			bb = square(b);
//			aa = square(a);
//
////			xyz[0] = r0 + SphericalGravityCenter[0];
////			index = getCellIndex(xyz);
////			rho = rhoField[index];
//			retcode = p->interpolateDensity(&rho, r); //g/cm**3
//			g = -SphericalGravityGetAt(r);
//			g_rho = g * rho;
//			P = (bb * (P1 + a * g_rho) - aa * (P2 + b * g_rho)) / (bb - aa);
////			TRACEF("%lld: P2 P1 P = %e %e %e ; r2 r1 r = %e %e %e ; g rho = %e %e", i, P2, P1, P, r2, r1, r, g, rho);
//
//			if(rr)
//			{
//				rr[i] = r;
//				gg[i] = g;
//				dd[i] = rho;
//			}
//
//			P_r[i] = P;
//			r2 = r1;
//			r1 = r;
//			P2 = P1;
//			P1 = P;
//		}
//
//		// P and rho have the outermost values from the P integration loop.
//		// Caclculate correction as if intergrating from the surface inwards.
//		const double K53 = 3.161128e+12;
//		const double K43 = 4.881986e+14;
//		double K = 3.828648955e+14;
//		P1 = K * POW(rho, Gamma);
//		P = P1 - P2; // add to P_r to get the right boundary conditions.
//		P = K * POW(2e9, Gamma);
////		TRACEF("Pressure shift: %e = %e - %e", P2, P1, P);
//		// Shift pressure and replace it with the gas energy.
//		for(size_t i = 0; i < SphericalGravityActualNumberOfBins; i++)
//		{
//			a = P_r[i] + P;
//			r = SphericalGravityBinLeftEdges[i];
//			retcode = p->interpolateDensity(&rho, r); //g/cm**3
//			b = a / rho / gammaMinusOne;
//			P_r[i] = a;
//			U_r[i] = b;
//		}
//
//		if(rr)
//		{
//			WriteInitialProfile("initialProfile", rr, dd, gg, P_r, U_r, SphericalGravityActualNumberOfBins, K,
//								(double)Gamma);
//		}
		for(int k = 0; k < GridDimension[2]; k++)
		{
			FLOAT z = CELLCENTER(2, k);
			FLOAT zz = square(z - SphericalGravityCenter[2]);
			for(int j = 0; j < GridDimension[1]; j++)
			{
				FLOAT y = CELLCENTER(1, j);
				FLOAT yy_zz = square(y - SphericalGravityCenter[1]) + zz;
				for(int i = 0; i <= GridDimension[0]; i++)
				{
					FLOAT x = CELLCENTER(0, i);
					FLOAT r = sqrt(square(x - SphericalGravityCenter[0]) + yy_zz);
					int index = ELT(i, j, k);

					float rho = rhoField[index];
					int retval;
					if(p->internalEnergyData)
					{
						retval = p->interpolateInternalEnergy(&U, r);
					}
					else
					{
						U = 0;
//						size_t rbin = SphericalGravityComputeBinIndex(r);
//						r1 = SphericalGravityBinLeftEdges[rbin];
//						r2 = SphericalGravityBinLeftEdges[rbin + 1];
//						// P_r[] has the specific gas energy now.
//						U1 = U_r[rbin];
//						U2 = U_r[rbin + 1];
//						U = ((r2 - r) * U1 + (r - r1) * U2) / (r2 - r1);
					}
					float gasE = U; // / rho / gammaMinusOne;
					if((1 + debug) && MyProcessorNumber == ROOT_PROCESSOR)
						if(j == (GridDimension[1] / 2) && k == (GridDimension[2] / 2) && i >= GridDimension[0] / 2)
							TRACEGF("i,j,k=%04d,%04d,%04d, x,y,z=(%4f,%4f,%4f), r[km]=%4f, rho=%e, U=%e", i, j, k,
									x * 1e-5, y * 1e-5, z * 1e-5, r * 1e-5, rho, gasE);

					totEField[index] = U; // add kinetic and magnetic energy below
				}
			}
		}
//
//		if(debug)
//		{
//			for(size_t i = 0; i < SphericalGravityActualNumberOfBins; i++)
//			{
//				r = SphericalGravityBinLeftEdges[i];
//				retcode = p->interpolateDensity(&rho, r); //g/cm**3
//				P = P_r[i];
//				U = U_r[i];
//				TRACEGF0("Intergrated ressure: %lld: r=%e P=%e U=%e rho=%e gamma-1=%f", i, r, P, U, rho, gammaMinusOne);
//			}
//		}
//
//		delete P_r;
//		delete U_r;
//		delete rr;
//		delete gg;
//		delete dd;
	}
	else
	{
		// If using gas energy field, the gas energy had been calculated
		// already, otherwise we need to calculate it here.
		if(hasGasEField)
		{
			arr_cpy(totEField, gasEField, gridSize);
		}
		else
		{
			for(float* TE = totEField; TE < totEField_end; TE++)
				*TE = EOSPolytropicFactor * POW(*RHO++, gammaMinusOne) / gammaMinusOne;
		}
	}

//Add kinetic energy
	for(float* TE = totEField; TE < totEField_end; TE++)
	{
		*TE += 0.5 * (square(*VX++) + square(*VY++) + square(*VZ++));
	}

//Add magnetic energy
	if(BX)
	{
		RHO = rhoField;
		float rho;
		for(float* TE = totEField; TE < totEField_end; TE++)
		{
			if((rho = *RHO++) > tiny_number)
				*TE += 0.5 * (square(*BX++) + square(*BY++) + square(*BZ++)) / rho;
			else
			{
				BX++;
				BY++;
				BZ++;
			}
		}
	}

// Boiler plate code:
//  if(DualEnergyFormalism )
//    for(index=0;index<size;index++)
//    BaryonField[ gesENum ][index] =
//       BaryonField[ totENum ][index]
//      - 0.5*(BaryonField[ vNum[0] ][index]*BaryonField[ vNum[0] ][index] +
//	     BaryonField[ vNum[1] ][index]*BaryonField[ vNum[1] ][index] +
//	     BaryonField[ vNum[2] ][index]*BaryonField[ vNum[2] ][index])
//      - 0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index] +
//             BaryonField[ByNum][index]*BaryonField[ByNum][index] +
//             BaryonField[BzNum][index]*BaryonField[BzNum][index])/BaryonField[ rhoNum ][index];

	if(HydroMethod == MHD_RK || usingVectorPotential)
	{
		UseMHDCT = FALSE;
		for(int dim = 0; dim < 3; dim++)
		{
			delete MagneticField[dim];
			delete ElectricField[dim];
			MagneticField[dim] = NULL;
			ElectricField[dim] = NULL;
		}
	}

	printf("Initialized grid, phase 2 (total energy.)\n");
	return SUCCESS;
}

int grid::WriteRadialProfile(char* name)
{
	if(MyProcessorNumber != ProcessorNumber)
		return FAIL;

	size_t num, ijk[3];
	FLOAT x, y, z, r;
	float rho, E, U, vx, vy, vz;
	float Bx = 0, By = 0, Bz = 0, rho_Ni = 0, g = 0, gamma = 0, P = 0;
	char filename[FILENAME_MAX];
	FLOAT xyz1[MAX_DIMENSION], xyz2[MAX_DIMENSION];
	FLOAT dx = CellWidth[0][0];
	size_t i1 = 0, i2 = 0;

	const int SAMPLE_OX = 0;
	const int SAMPLE_DIAG = 1;
	int sampleWhere = SAMPLE_DIAG;
	switch(sampleWhere)
	{
	case SAMPLE_OX:
		// Sample a row of cells along x,
		// (0,0,0)..(DomainRightEdge[0], 0, 0)
		arr_set(xyz1, MAX_DIMENSION, dx * 1e-6);
		arr_set(xyz2, MAX_DIMENSION, dx * 1e-6);
		xyz2[0] = DomainRightEdge[0] - dx * 1e-6;
		break;
	case SAMPLE_DIAG:
	default:
		// Sample cells from the center along the diagonal,
		// (0,0,0)..(DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2]):
		arr_set(xyz1, MAX_DIMENSION, dx * 1e-6);
		arr_set(xyz2, MAX_DIMENSION, -dx * 1e-6);
		arr_xpy(xyz2, DomainRightEdge, MAX_DIMENSION);
	}

//Check if the segment xyz1..xyz2 is inside the grid.
	for(int dim = 0; dim < GridRank; dim++)
	{
		if(xyz1[dim] < GridLeftEdge[dim] || (GridRightEdge[dim] < xyz2[dim]))
			return FAIL;
	}

	switch(sampleWhere)
	{
	case SAMPLE_OX:
	case SAMPLE_DIAG:
	default:
		get_ijk_index(ijk, xyz1);
		i1 = ijk[0];
		get_ijk_index(ijk, xyz2);
		i2 = ijk[0];
	}

	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	int totENum, rhoNum, vxNum, vyNum, vzNum;
	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	float *totEField, *rhoField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	if(IdentifyPhysicalQuantities(rhoNum, gasENum, vxNum, vyNum, vzNum, totENum, BxNum, ByNum, BzNum) == FAIL)
	{
		ENZO_FAIL("grid::WriteRadialProfile: Error in IdentifyPhysicalQuantities for UseMHD==true.");
	}

	if(UseMHD)
	{
		BxField = BaryonField[BxNum];
		ByField = BaryonField[ByNum];
		BzField = BaryonField[BzNum];
	}

	if(UseBurning)
	{
		int rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
		if(rhoNiNum < 0)
			ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
		rhoNiField = BaryonField[rhoNiNum];
	}

	rhoField = BaryonField[rhoNum];
	vxField = BaryonField[vxNum];
	vyField = BaryonField[vyNum];
	vzField = BaryonField[vzNum];
	totEField = BaryonField[totENum];
	if(hasGasEField)
		gasEField = BaryonField[gasENum];

	for(int format = 0; format < 2; format++)
	{
		switch(format)
		{
		case 0:
			sprintf(filename, "%s.radialProfile", name);

			break;
		case 1:
			sprintf(filename, "%s.radialProfile.py", name);
			break;
		}

		FILE* file = fopen(filename, "w");
		if(file == NULL)
		{
			ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
		}

		switch(format)
		{
		case 0:
			fprintf(file, "%f seconds\n", Time);
			break;
		case 1:
			break;
		}

		switch(format)
		{
		case 0:
			fprintf(file, "# rownum   x y z   rho rho_Ni   E U   vx vy vz   Bx By Bz   g   gamma P\n");
			break;
		case 1:
			fprintf(file, "type('', (), dict(\n");
			fprintf(file, "time = %f,\n", Time);
			fprintf(file, "cols = 'rownum x y z rho rho_Ni E U vx vy vz Bx By Bz g gamma P'.split(),\n");
//			fprintf(file, "radialProfile.dtype = np.dtype([\n"
//				"  ('rownum',  int),\n"
//				"  ('x'    , float),\n"
//				"  ('y'    , float),\n"
//				"  ('z'    , float),\n"
//				"  ('rho'  , float),\n"
//				"  ('rhoNi', float),\n"
//				"  ('Ni'   , float),\n"
//				"  ('E'    , float),\n"
//				"  ('U'    , float),\n"
//				"  ('vx'   , float),\n"
//				"  ('vy'   , float),\n"
//				"  ('vz'   , float),\n"
//				"  ('Bx'   , float),\n"
//				"  ('By'   , float),\n"
//				"  ('Bz'   , float),\n"
//				"  ('g'    , float),\n"
//				"  ('gamma', float),\n"
//				"  ('P'    , float),\n"
//				"])\n");
			fprintf(file, "data = np.array([\n");
			break;
		}

		arr_set(ijk, MAX_DIMENSION, 0);
		for(size_t i = i1; i <= i2; i++)
		{
			num = i - i1 + 1;

			switch(sampleWhere)
			{
			case SAMPLE_OX:
				ijk[0] = i;
				break;
			case SAMPLE_DIAG:
			default:
				arr_set(ijk, GridRank, i);
			}

			size_t index = getCellIndex(ijk);

			get_xyz(xyz1, ijk);
			x = xyz1[0];
			y = xyz1[1];
			z = xyz1[2];
			r = distancel(xyz1, SphericalGravityCenter, GridRank);
			E = totEField[index];
			rho = rhoField[index];
			vx = vxField[index];
			vy = vyField[index];
			vz = vzField[index];
			gamma = Gamma;

			if(BxField)
			{
				Bx = BxField[index];
				By = ByField[index];
				Bz = BzField[index];
			}

			if(gasEField)
				U = gasEField[index];
			else
			{
				U = square(Bx) + square(By) + square(Bz);
				U /= rho;
				U += square(vx) + square(vy) + square(vz);
				U *= -0.5;
				U += E;
			}
			P = (gamma - 1) * U * rho;

			if(UseBurning)
				rho_Ni = rhoNiField[index];

			if(UseSphericalGravity)
			{
				g = SphericalGravityGetAt(r);
			}

//		fprintf(file, "num   x   y   z   rho   E   U   vx   vy   vz   Bx   By   Bz   g   gamma P\n");
			switch(format)
			{
			case 0:
				fprintf(file, "%lld   %e %e %e   %e %e   %e %e   %e %e %e   %e %e %e   %e   %f %e\n", /**/
						num, x, y, z, rho, rho_Ni, E, U, vx, vy, vz, Bx, By, Bz, g, gamma, P);
				break;
			case 1:
				fprintf(file, "[ %lld ,   %e , %e , %e ,   %e , %e ,   %e , %e ,   %e , %e , %e ,   %e , %e , %e"
						" ,   %e ,   %f , %e ],\n",
						num, x, y, z, rho, rho_Ni, E, U, vx, vy, vz, Bx, By, Bz, g, gamma, P);
				break;
			}
		}

		switch(format)
		{
		case 0:
			break;
		case 1:
			fprintf(file, "])))\n");
			break;
		} // end for data

		if(file != stderr && file != stdout)
		{
			fclose(file);
			fprintf(stderr, "'%s' written.\n", filename);
		}
	} // end for format

	return SUCCESS;
}
