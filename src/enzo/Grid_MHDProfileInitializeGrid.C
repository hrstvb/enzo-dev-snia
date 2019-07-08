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
#include <stdlib.h>
#include <cmath>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "mylimiters.h"
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
#include "TriSpherePerturb.C"

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

inline float internalEnergy(float rho, float rhoNi, MHDInitialProfile* profile, FLOAT radius)
{
	double T = 0, E = 0;
	if(profile->interpolateTemperature(&T, radius))
		T = (rhoNi / rho > 0.5) ? 5e9 : 1e6;
	E += 1.5 * (1 - 0.75 * rhoNi / rho) * R_gas * T / 14;
	E += EOSPolytropicFactor * pow(rho, Gamma - 1) / (Gamma - 1);
	return E;
}

inline float internalEnergy(float* rhoField, float* rhoNiField, long long index, MHDInitialProfile* profile,
FLOAT radius)
{
	return internalEnergy(rhoField[index], rhoNiField[index], profile, radius);
}

int grid::MHD_SNIA_GetFields(float** densityField, float** totalEnergyField, float** internalEnergyField,
float** vxField, float** vyField, float** vzField, float** vFields,
float** BxField, float** ByField, float** BzField, float** BFields,
float** burnedDensityField,
float** phiField,
float** gravPotentialField)
{
	int rhoNum, gasENum, vxNum, vyNum, vzNum, totENum;
	int BxNum, ByNum, BzNum, phiNum;
	int rhoNiNum;
	int gravPotentialNum;

	// IdentifyPhysicalQuantities set gasENum=0 if not having Gas Energy field.
	bool hasGasEField = DualEnergyFormalism; // && (HydroMethod != Zeus_Hydro); //[BH]
	BxNum = ByNum = BzNum = phiNum = -1; // IdentifyPhysicalQuantities doesn't change these args if (not UseMHD).
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
	else
	{
		if(BxField)
			*BxField = NULL;
		if(ByField)
			*ByField = NULL;
		if(BzField)
			*BzField = NULL;
		if(BFields)
			BFields[0] = BFields[1] = BFields[2] = NULL;
	}

	if(burnedDensityField)
	{
		if(UseBurning)
		{
			rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
			if(rhoNiNum < 0)
				ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
			*burnedDensityField = BaryonField[rhoNiNum];
		}
		else
		{
			*burnedDensityField = NULL;
		}
	}

	if(phiField)
		*phiField = (phiNum >= 0) ? BaryonField[phiNum] : NULL;

	if(gravPotentialField)
	{
		if(WritePotential)
		{
			gravPotentialNum = FindField(GravPotential, FieldType, NumberOfBaryonFields);
			if(gravPotentialNum < 0)
				ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for GravPotential.")
			*gravPotentialField = BaryonField[gravPotentialNum];
		}
		else
		{
			*gravPotentialField = NULL;
		}
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
		FieldType[NumberOfBaryonFields++] = Density_56Ni; //[BH]
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

	TRACEGF("INITIALIZING GRID PHASE 1, FIELDS.");

	this->AllocateGrids();

//	int totENum, rhoNum, vxNum, vyNum, vzNum;
//	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField,
						NULL,
						&rhoNiField, NULL, NULL);

	int debugnanflag = 0; //[BH]
	float gammaMinusOne = Gamma - 1;

	// Initialize density and a sphere of Nickel density.
	for(int k = 0; k < GridDimension[2]; k++)
	{
		FLOAT z = CELLCENTER(2, k);
		FLOAT rz = z - SphericalGravityCenter[2];
		FLOAT zz = square(rz);
		for(int j = 0; j < GridDimension[1]; j++)
		{
			FLOAT y = CELLCENTER(1, j);
			FLOAT ry = y - SphericalGravityCenter[1];
			FLOAT yy_zz = square(ry) + zz;
			for(int i = 0; i <= GridDimension[0]; i++)
			{
				int index = ELT(i, j, k);
				float vr, vx, vy, vz, vv, Bx, By, Bz, BB;

				FLOAT x = CELLCENTER(0, i);
				FLOAT rx = x - SphericalGravityCenter[0];
				FLOAT r = sqrt(square(rx) + yy_zz);

				float rho, rhoNi = 0, T = 0;
				int retcode = p->interpolateDensity(&rho, r); //g/cm**3
				//retcode = profileInterpolate(&T, temperatureData, r, radiusData, p.nRows, radiusSO); //K
				bool isBurned = false;
				double r1, r2;

				switch(PerturbationMethod)
				{
				case 0:
					isBurned = r <= InitialBurnedRadius;
					break;
				case 1: // Burned sphere, smooth perturb zone.
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
					{
						rhoNi = (r - r1) / PerturbationAmplitude;
						rhoNi *= 4 / (1 + 3 * rhoNi);
						rhoNi *= rho;
					}
					break;
				case 2: // Burned sphere, randomly burned perturb zone
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
						isBurned = rand() % 2;
					break;
				case 3: // Burned sphere, randomly burned perturb zone
						// Probability decreasing with radius.
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
					{
						isBurned = drand48() >= (r - r1) / PerturbationAmplitude;
					}
					break;
				case 4:
					isBurned = r <= InitialBurnedRadius;
//					if(PertrubationTopDensity > 0 && InitialBurnedRadius < r
//							&& r < InitialBurnedRadius + PerturbationAmplitude)
//						rho = PertrubationTopDensity;
					break;
				}

				if(isBurned)
					rhoNi = rho;
				rhoField[index] = rho;
				if(UseBurning)
					rhoNiField[index] = rhoNi;

//					if(debug + 1 && MyProcessorNumber == ROOT_PROCESSOR)
//						if((j == 0 || j == (GridDimension[1] - 1)) && (k == 0 || k == (GridDimension[2] - 1))
//								&& fabs(y) < fabs(GridRightEdge[1] - GridLeftEdge[0]) / 4
//								&& fabs(z) < fabs(GridRightEdge[2] - GridLeftEdge[2]) / 4)
////							if(j == (GridDimension[1]) && k == (GridDimension[2]) && i >= GridDimension[0])
//							TRACEF("i,j,k=%03d,%04d,%04d, x,y,z=(%4f,%4f,%4f), rx,ry,rz=(%4f,%4f,%4f), r=%4f, " //
//							"rho=%e, U=%e, E=%e, v_r=%e, v^2=%e, B^2=%e, "//
//							"burned=%d, K=%e, gamma-1=%f",//
//									i, j, k, x * 1e-5, y * 1e-5, z * 1e-5, rx * 1e-5, ry * 1e-5, rz * 1e-5, r * 1e-5, //
//									rho, gasE, totE, vr, vv, BB, //
//									isBurned, EOSPolytropicFactor, gammaMinusOne);
			} //end baryonfield initialize
		}
	}

	if(PerturbationMethod == 4 && burnedRadius > 0)
	{
		FILE* fptr = NULL;
		if(MyProcessorNumber == ROOT_PROCESSOR && isTopGrid())
			fptr = fopen("trisphere-data.py", "w");

		if(fptr)
		{
			fprintf(fptr, "type('', (), dict(\n");

			fprintf(fptr, ""
					"## Unperturbed sphere (primary surface) triangulation parameters:\n"
					"R = %e , # Primary radius\n"
					"nV = %lld , # Number of vertices\n"
					"nE = %lld , # Number of edges\n"
					"nF = %lld , # Number of faucets \n"
					"## Perturbations parameters:\n"
					"A = %e , # Approx. perturbation height\n"
					"bottomBaseSize = %4.2f , # relative to the primary facet (linear)\n"
					"topBaseSize = %4.2f , # relative to the primary facet (linear)\n",
					triSphere->R, triSphere->nV, triSphere->nE, triSphere->nF, triSphere->A, triSphere->bottomBaseSize,
					triSphere->topBaseSize);
			fprintf(fptr, "## Grid parameters:\n"
					"gridLevel = %lld , # Hierarchy level\n"
					"gridID = %lld , # grid ID in level\n"
					"grid_dx = %e , # Grid spatial step\n"
					"numGhostZones = %lld , \n"
					"gridEdges = [\n"
					"  [ %e , %e , %e ], # grid left edges\n"
					"  [ %e , %e , %e ]], #grid right edges\n",
					0, ID, TopGridDx[0], NumberOfGhostZones, GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2],
					GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);

			fprintf(fptr, "\n# Perturbation data for the grid:\ndata=[\n\n");
		}

//Perturb in the shell between r1 < r <r2.
// Calculate position with respect to each side.
		for(size_t k_face = 0; k_face < triSphere->nF; k_face++)
//		for(size_t k_face = 9; k_face < 10; k_face++)
		{
			long lijk[3], rijk[3];
			double *ledge, *redge;
			FLOAT xyz[3], xyz2[3];
			size_t numInside = 0;

			triSphere->getSpikeEnclosingRectangle(k_face, &ledge, &redge);
			if(fptr)
			{
				fprintf(fptr, "[ %lld ,\n", k_face);
				triSphere->fprint_cache(fptr);
				fprintf(fptr, "],\n[");
			}

			if(0 > intersect(ledge, redge))
			{
				if(fptr)
					fprintf(fptr, "]],\n");
				continue;
			}

			get_ijk_index(lijk, ledge);
			get_ijk_index(rijk, redge);

			for(size_t k = lijk[2]; k < rijk[2]; k++)
			{
				xyz[2] = CELLCENTER(2, k);
				//xyz2[2] = CELLCENTER(2, k);
				for(size_t j = lijk[1]; j < rijk[1]; j++)
				{
					xyz[1] = CELLCENTER(1, j);
					//xyz2[1] = CELLCENTER(1, j);
					bool wasInside = false;
					for(size_t i = lijk[0]; i < rijk[0]; i++)
					{
						xyz[0] = CELLCENTER(0, i);
						//xyz2[0] = CELLCENTER(0, i);

						bool isInsideSpike = triSphere->isInSpike(k_face, xyz);
						if(!isInsideSpike)
							continue;

						if(fptr)
							fprintf(fptr, "[ %e , %e , %e ],\n", xyz[0], xyz[1], xyz[2]);

						size_t index = ELT(i, j, k);
						double DX = TopGridDx[0] / CellWidth[0][0];
						double massfrac = (DX > .75) ? 1 : (DX > .4) ? .2 : .01;
						massfrac = 4 * massfrac / (1 + 3 * massfrac);
						massfrac = 1;
						float rho;
						if(PertrubationBottomDensity > 0)
							rhoField[index] = rho = PertrubationBottomDensity;
						else
							rho = rhoField[index];
						rhoNiField[index] = massfrac * rho;
						wasInside = true;
						numInside++;
					}
				}
			} // end i,j,k loops

			if(fptr)
				fprintf(fptr, "]],\n");
		} // end k_face loop

		if(fptr)
		{
			fprintf(fptr, "\n] #end data\n");
			fprintf(fptr, ") #end dict\n");
			fprintf(fptr, ") #end type\n\n");
			fclose(fptr);
			fptr == NULL;
		}
	}

	//Initialize all other fields.
	for(int k = 0; k < GridDimension[2]; k++)
	{
		FLOAT z = CELLCENTER(2, k);
		FLOAT rz = z - SphericalGravityCenter[2];
		FLOAT zz = square(rz);
		for(int j = 0; j < GridDimension[1]; j++)
		{
			FLOAT y = CELLCENTER(1, j);
			FLOAT ry = y - SphericalGravityCenter[1];
			FLOAT yy_zz = square(ry) + zz;
			for(int i = 0; i <= GridDimension[0]; i++)
			{
				int index = ELT(i, j, k);
				float vr, vx, vy, vz, vv, Bx, By, Bz, BB;

				//!strcmp(p->profileType, "RADIAL"))
				{
					FLOAT x = CELLCENTER(0, i);
					FLOAT rx = x - SphericalGravityCenter[0];
					FLOAT r = sqrt(square(rx) + yy_zz);

					float rho = rhoField[index];
					float rhoNi = rhoNiField[index];
					float T = 0;

					vx = GridVelocity[0];
					vy = GridVelocity[1];
					vz = GridVelocity[2];
					if(p->radialVelocityData)
					{
						p->interpolateRadialVelocity(&vr, r);
						vx += vr * rx / r;
						vy += vr * ry / r;
						vz += vr * rz / r;
					}
					vx = vy = vz = 0;
//					vz=1e5;
//					if((i==30 || i==29) && j == 10 && k == 11)
//						vy = 1000e5;

//					{
//						for(int i = 0; i<p->nRowsInternalEnergy; i++)
//							printf(" %e, ", p->internalEnergyData[i]);
//						printf("\n");
//					}
					float gasE = 0;
					if(p->internalEnergyData)
					{
						p->interpolateInternalEnergy(&gasE, r); //g/cm**3
					}
					else
					{
//						gasE = EOSPolytropicFactor * POW(rho, gammaMinusOne) / gammaMinusOne;
						gasE = internalEnergy(rho, rhoNi, p, r);
					}

//					if(UseBurning && isBurned && burnedRadiusPE)
//					{
//						rho *= gasE / (gasE + BurningEnergyRelease);
//						gasE += BurningEnergyRelease;
//					}

					float totE = gasE;
//					vv = square(vx) * square(vy) + square(vz);
					vv = vx * vx + vy * vy + vz * vz;
					totE += 0.5 * vv;
					if(BxField)
					{
						MHDProfileInitExactB(&Bx, &By, &Bz, x, y, z);
						BB = square(Bx) + square(By) + square(Bz);
						totE += 0.5 * BB / rho;
						BxField[index] = Bx;
						ByField[index] = By;
						BzField[index] = Bz;
					}

//					rhoField[index] = rho;
					if(hasGasEField)
						gasEField[index] = gasE;
					totEField[index] = totE;
					vxField[index] = vx;
					vyField[index] = vy;
					vzField[index] = vz;
//					if(UseBurning)
//						rhoNiField[index] = rhoNi;

//						if(debug + 1 && MyProcessorNumber == ROOT_PROCESSOR)
//							if((j == 0 || j == (GridDimension[1] - 1)) && (k == 0 || k == (GridDimension[2] - 1))
//									&& fabs(y) < fabs(GridRightEdge[1] - GridLeftEdge[0]) / 4
//									&& fabs(z) < fabs(GridRightEdge[2] - GridLeftEdge[2]) / 4)
////							if(j == (GridDimension[1]) && k == (GridDimension[2]) && i >= GridDimension[0])
//								TRACEF("i,j,k=%03d,%04d,%04d, x,y,z=(%4f,%4f,%4f), rx,ry,rz=(%4f,%4f,%4f), r=%4f, " //
//								"rho=%e, U=%e, E=%e, v_r=%e, v^2=%e, B^2=%e, "//
//								"burned=%d, K=%e, gamma-1=%f",//
//										i, j, k, x * 1e-5, y * 1e-5, z * 1e-5, rx * 1e-5, ry * 1e-5, rz * 1e-5,
//										r * 1e-5, //
//										rho, gasE, totE, vr, vv, BB, //
//										isBurned, EOSPolytropicFactor, gammaMinusOne);
				}
			} //end baryonfield initialize
		}
	}

//	for(int k = 0; k < GridDimension[2]; k++)
//	{
//		FLOAT z = CELLCENTER(2, k);
//		FLOAT rz = z - SphericalGravityCenter[2];
//		FLOAT zz = square(rz);
//		for(int j = 0; j < GridDimension[1]; j++)
//		{
//			FLOAT y = CELLCENTER(1, j);
//			FLOAT ry = y - SphericalGravityCenter[1];
//			FLOAT yy_zz = square(ry) + zz;
//			for(int i = 0; i <= GridDimension[0]; i++)
//			{
//				int index = ELT(i, j, k);
//				float vr, vx, vy, vz, vv, Bx, By, Bz, BB;
//
//				//!strcmp(p->profileType, "RADIAL"))
//				{
//					FLOAT x = CELLCENTER(0, i);
//					FLOAT rx = x - SphericalGravityCenter[0];
//					FLOAT r = sqrt(square(rx) + yy_zz);
//
//					float rho = rhoField[index];
//					float rhoNi = rhoNiField[index];
//					float T = 0;
//
//					vx = GridVelocity[0];
//					vy = GridVelocity[1];
//					vz = GridVelocity[2];
//					if(p->radialVelocityData)
//					{
//						p->interpolateRadialVelocity(&vr, r);
//						vx += vr * rx / r;
//						vy += vr * ry / r;
//						vz += vr * rz / r;
//					}
//					vx = vy = vz = 0;
////					vz=1e5;
////					if((i==30 || i==29) && j == 10 && k == 11)
////						vy = 1000e5;
//
////					{
////						for(int i = 0; i<p->nRowsInternalEnergy; i++)
////							printf(" %e, ", p->internalEnergyData[i]);
////						printf("\n");
////					}
//					float gasE = 0;
//					if(p->internalEnergyData)
//					{
//						p->interpolateInternalEnergy(&gasE, r); //g/cm**3
//					}
//					else
//					{
////						gasE = EOSPolytropicFactor * POW(rho, gammaMinusOne) / gammaMinusOne;
//						gasE = internalEnergy(rho, rhoNi, p, r);
//					}
//
////					if(UseBurning && isBurned && burnedRadiusPE)
////					{
////						rho *= gasE / (gasE + BurningEnergyRelease);
////						gasE += BurningEnergyRelease;
////					}
//
//					float totE = gasE;
////					vv = square(vx) * square(vy) + square(vz);
//					vv = vx * vx + vy * vy + vz * vz;
//					totE += 0.5 * vv;
//					if(BxField)
//					{
//						MHDProfileInitExactB(&Bx, &By, &Bz, x, y, z);
//						BB = square(Bx) + square(By) + square(Bz);
//						totE += 0.5 * BB / rho;
//						BxField[index] = Bx;
//						ByField[index] = By;
//						BzField[index] = Bz;
//					}
//
////					rhoField[index] = rho;
//					if(hasGasEField)
//						gasEField[index] = gasE;
//					totEField[index] = totE;
//					vxField[index] = vx;
//					vyField[index] = vy;
//					vzField[index] = vz;
////					if(UseBurning)
////						rhoNiField[index] = rhoNi;
//
////						if(debug + 1 && MyProcessorNumber == ROOT_PROCESSOR)
////							if((j == 0 || j == (GridDimension[1] - 1)) && (k == 0 || k == (GridDimension[2] - 1))
////									&& fabs(y) < fabs(GridRightEdge[1] - GridLeftEdge[0]) / 4
////									&& fabs(z) < fabs(GridRightEdge[2] - GridLeftEdge[2]) / 4)
//////							if(j == (GridDimension[1]) && k == (GridDimension[2]) && i >= GridDimension[0])
////								TRACEF("i,j,k=%03d,%04d,%04d, x,y,z=(%4f,%4f,%4f), rx,ry,rz=(%4f,%4f,%4f), r=%4f, " //
////								"rho=%e, U=%e, E=%e, v_r=%e, v^2=%e, B^2=%e, "//
////								"burned=%d, K=%e, gamma-1=%f",//
////										i, j, k, x * 1e-5, y * 1e-5, z * 1e-5, rx * 1e-5, ry * 1e-5, rz * 1e-5,
////										r * 1e-5, //
////										rho, gasE, totE, vr, vv, BB, //
////										isBurned, EOSPolytropicFactor, gammaMinusOne);
//				}
//			} //end baryonfield initialize
//		}
//	}

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
	TRACEGF("... DONE INITIALIZING GRID FULL PHASE 1.");

	return SUCCESS;
}

/*
 * Keeps the burned fraction equal to 1 for the initial burned raadius.
 */
int grid::MHDSustainInitialBurnedRegionGrid()
{
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
	FLOAT lxyz[MAX_DIMENSION], rxyz[MAX_DIMENSION];
	long lindex[MAX_DIMENSION], rindex[MAX_DIMENSION];

	arr_set(lxyz, MAX_DIMENSION, 0);
	arr_set(rxyz, MAX_DIMENSION, 0);
	for(int dim = 0; dim < GridRank; dim++)
	{
		lxyz[dim] = SphericalGravityCenter[dim] - InitialBurnedRadius;
		rxyz[dim] = SphericalGravityCenter[dim] + InitialBurnedRadius;
	}

	if(intersect(lxyz, rxyz))
		return SUCCESS;

	get_ijk_index(lindex, lxyz);
	get_ijk_index(rindex, rxyz);

	const FLOAT RR = square(InitialBurnedRadius);
	const FLOAT X0 = SphericalGravityCenter[0];
	switch(GridRank)
	{
	case 1:
	{
		FLOAT xmin = X0 - RR;
		FLOAT xmax = X0 + RR;
		for(int i = lindex[0]; i <= rindex[0]; i++)
		{
			FLOAT x = CELLCENTER(0, i);
			if(xmin <= x && x <= xmax)
				rhoNiField[i] = rhoField[i];
		}
		break;
	}
	case 2:
		for(int j = lindex[1]; j < rindex[1]; j++)
		{
			FLOAT xxmax = RR - square(CELLCENTER(1, j) - SphericalGravityCenter[1]);
			if(xxmax < 0)
				continue;
			size_t index = ELT(lindex[0], j, 0);
			for(int i = lindex[0]; i <= rindex[0]; i++)
			{
				FLOAT xx = square(CELLCENTER(0, i) - X0);
				if(xx <= xxmax)
					rhoNiField[index] = rhoField[index];
				index++;
			}
		}
		break;

	case 3:
		for(int k = lindex[2]; k <= rindex[2]; k++)
		{
			FLOAT zz = square(CELLCENTER(2, k) - SphericalGravityCenter[2]);
			for(int j = lindex[1]; j < rindex[1]; j++)
			{
				FLOAT xxmax = RR - square(CELLCENTER(1, j) - SphericalGravityCenter[1]) - zz;
				if(xxmax < 0)
					continue;
				size_t index = ELT(lindex[0], j, k);
				for(int i = lindex[0]; i <= rindex[0]; i++)
				{
					FLOAT xx = square(CELLCENTER(0, i) - X0);
					if(xx <= xxmax)
						rhoNiField[index] = rhoField[index];
					// Debug code:
					// rhoNiField[index] = (xx <= xxmax) ? 0 : rhoField[index];
					index++;
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
						&BzField,
						NULL,
						&rhoNiField, NULL, NULL);

	float gammaMinusOne = Gamma - 1;
	size_t gridSize = GetGridSize();
	const float* totEField_end = totEField + gridSize;

	float* RHO = rhoField;
	float* NI = rhoNiField;
	float* GE = gasEField;
	float* VX = vxField;
	float* VY = vyField;
	float* VZ = vzField;
	float* BX = BxField;
	float* BY = ByField;
	float* BZ = BzField;

	FLOAT z, zz, y, yy_zz, x, r;
	float gasE;
	size_t index;
	if(hasGasEField)
	{
// If using gas energy field, the gas energy had been calculated
// already, otherwise we need to calculate it here.
		arr_cpy(totEField, gasEField, gridSize);
	}
	else
	{
		bool hasInternalEnergyProfile = InitRadialPressureFromCentral && p->internalEnergyData;
		for(int k = 0; k < GridDimension[2]; k++)
		{
			z = CELLCENTER(2, k);
			zz = square(z - SphericalGravityCenter[2]);
			for(int j = 0; j < GridDimension[1]; j++)
			{
				y = CELLCENTER(1, j);
				yy_zz = square(y - SphericalGravityCenter[1]) + zz;
				index = ELT(0, j, k);

				for(int i = 0; i <= GridDimension[0]; i++)
				{
					x = CELLCENTER(0, i);
					r = sqrt(square(x - SphericalGravityCenter[0]) + yy_zz);

					if(hasInternalEnergyProfile)
						p->interpolateInternalEnergy(&gasE, r);
					else
						gasE = internalEnergy(RHO, NI, index, p, r);

					totEField[index] = gasE;

					if((1 + debug) && MyProcessorNumber == ROOT_PROCESSOR)
						if(j == (GridDimension[1] / 2) && k == (GridDimension[2] / 2) && i >= GridDimension[0] / 2)
							TRACEGF("i,j,k=%04d,%04d,%04d, x,y,z=(%4f,%4f,%4f), r[km]=%4f, rho=%e, U=%e", i, j, k,
									x * 1e-5, y * 1e-5, z * 1e-5, r * 1e-5, rho, gasE);

					index++;
				}
			}
		}
	}

	//special pressure init test for monogrid.
	if(1)
	{
		TRACE
		;

#define DEBUG_HSE 0
#define SET_HSE(i,j,k,k2,E) do{ \
/*		TRACEF("%lld %lld %lld", (00)   + (i), (k2)-1 - (j), (k2)-1 - (k)); */ \
/*		TRACEF("%lld %lld %lld", (k2)-1 - (i), (k2)-1 - (j), (k2)-1 - (k)); */ \
/*		TRACEF("%lld %lld %lld", (00)   + (i), (00)   + (j), (k2)-1 - (k)); */ \
/*		TRACEF("%lld %lld %lld", (k2)-1 - (i), (00)   + (j), (k2)-1 - (k)); */ \
/*		TRACEF("%lld %lld %lld", (00)   + (i), (k2)-1 - (j), (00)   + (k)); */ \
/*		TRACEF("%lld %lld %lld", (k2)-1 - (i), (k2)-1 - (j), (00)   + (k)); */ \
/*		TRACEF("%lld %lld %lld", (00)   + (i), (00)   + (j), (00)   + (k)); */ \
/*		TRACEF("%lld %lld %lld", (k2)-1 - (i), (00)   + (j), (00)   + (k)); */ \
								totEField[ELT((00)   + (i), (k2)-1 - (j), (k2)-1 - (k))] += (E); \
								totEField[ELT((k2)-1 - (i), (k2)-1 - (j), (k2)-1 - (k))] += (E); \
								totEField[ELT((00)   + (i), (00)   + (j), (k2)-1 - (k))] += (E); \
								totEField[ELT((k2)-1 - (i), (00)   + (j), (k2)-1 - (k))] += (E); \
								totEField[ELT((00)   + (i), (k2)-1 - (j), (00)   + (k))] += (E); \
								totEField[ELT((k2)-1 - (i), (k2)-1 - (j), (00)   + (k))] += (E); \
								totEField[ELT((00)   + (i), (00)   + (j), (00)   + (k))] += (E); \
								totEField[ELT((k2)-1 - (i), (00)   + (j), (00)   + (k))] += (E); \
							}while(0)
//
		const size_t FIELD_SIZE = GetGridSize();
		const size_t K2 = GridDimension[2];
		const size_t K0 = K2 / 2;
		const FLOAT dx = CellWidth[0][0];
		FLOAT x, y, z;
		double P_c, P0, P1, P2, dP01, dP12;
		double r, rr, rho, rho_c, xi, xifactor;
		int err;

		arr_set(totEField, FIELD_SIZE, 0);

//Initialize the center cells
		double KPolytropic = EOSPolytropicFactor;
		double nPolytropic = 1 / (Gamma - 1);
		p->interpolateDensity(&rho_c, 0);
		P_c = KPolytropic * POW(rho_c, Gamma);
		TRACEF("%e %e %e %f", P_c, KPolytropic, rho_c, Gamma);
		xifactor = (nPolytropic + 1) * KPolytropic * POW(rho_c, 1 / nPolytropic - 1);
		xifactor /= 4 * M_PI * SphericalGravityConstant;
		xifactor = sqrt(xifactor);

		TRACEF("%lld .. %lld", K0, K2);

		for(size_t k = K0; k < K0 + 2; k++)
		{
			z = CELLCENTER(2, k) - SphericalGravityCenter[2];
			for(size_t j = K0; j < K0 + 2; j++)
			{
				y = CELLCENTER(1, j) - SphericalGravityCenter[1];
				for(size_t i = K0; i < K0 + 2; i++)
				{
					x = CELLCENTER(0, i) - SphericalGravityCenter[0];
					r = lenl(x, y, z);
					xi = r / xifactor;
					P0 = P_c * (1 - (nPolytropic + 1) / 6 * xi * xi);
					if(DEBUG_HSE)
						P0 = 1e1;
					TRACEF("ijk=%lld %lld %lld xyz=%e %e %e r,xi,r/xi=%e %e %e   P=%e", i, j, k, x, y, z, r, xi,
							xifactor, P0);
					SET_HSE(i, j, k, K2, P0);
				}
			}
		}

		if(1)
			// Integrate from the center cells out along z.
			for(size_t j = K0; j < K0 + 2; j++)
			{
				y = CELLCENTER(1, j) - SphericalGravityCenter[1];
				for(size_t i = K0; i < K0 + 2; i++)
				{
					x = CELLCENTER(0, i) - SphericalGravityCenter[0];
					P0 = totEField[ELT(i, j, K0)];
					P1 = totEField[ELT(i, j, K0 + 1)];
					dP01 = P1 - P0;
					for(size_t k = K0 + 2; k < K2; k++)
					{
						z = CELLCENTER(2, k-1) - SphericalGravityCenter[2];
						r = lenl(x, y, z);
						p->interpolateDensity(&rho, r);
						double g = -z / r * SphericalGravityGetAt(r);
						double slope = g * rho * dx;
						dP12 = invLimiter(dP01, slope, &err);
						P2 = P1 + dP12;
//						P0 = 2 * (dP01 * dP12 * dP12 * dP12 * dP12 + dP01 * dP01 * dP01 * dP01 * dP12)
//								/ ((dP01 * dP01 + dP12 * dP12) * (dP01 * dP01 + dP12 * dP12));
//					TRACEF("dP01=%e dP12=%e slope_in=%e slope_out=%e", dP01, dP12, slope, P0);
//					if(err)
						if(i == K0 && (j == K0 || j == K0 + 1))
						{
							TRACEF("i,j,k= %lld %lld %lld   K0,K2= %lld %lld    P0,P1,dP01= %e %e %e    slope,slope/dP01= %e %e    dP12,P2= %e %e",
									i, j, k, K0, K2, P0, P1, dP01, slope, slope / dP01, dP12, P2);
						}
						if(DEBUG_HSE)
							P2 = 1e2;
						SET_HSE(i, j, k, K2, P2);

						dP01 = dP12;
						P0 = P1;
						P1 = P2;
					}
				}
			}
		TRACEF("P_boundary=%e", P2);

		if(1)
			// Integrate from the center cells out along y.
			for(size_t k = K0; k < K2; k++)
			{
				z = CELLCENTER(2, k) - SphericalGravityCenter[2];
				for(size_t i = K0; i < K0 + 2; i++)
				{
					x = CELLCENTER(0, i) - SphericalGravityCenter[0];
					P0 = totEField[ELT(i, K0, k)];
					P1 = totEField[ELT(i, K0 + 1, k)];
					dP01 = P1 - P0;
//					dP01 = -abs(dP01);
//					if(i == K0 && k == K0 + 2)
//						TRACEF("P1 P0 dP01 %e %e %e", P1, P0, dP01);
					for(size_t j = K0 + 2; j < K2; j++)
					{
						y = CELLCENTER(1, j-1) - SphericalGravityCenter[1];
						r = lenl(x, y, z);
						p->interpolateDensity(&rho, r);
						double g = -y / r * SphericalGravityGetAt(r);
						double slope = g * rho * dx;
						err = 0 && (i == 56 && j == 109 && k == 57);
						if(err)
						{
//							TRACEF("K0, K2, i, j, k = %lld, %lld, %lld, %lld, %lld", K0, K2, i, j, k);
//							TRACEF("ijk=%lld %lld %lld, xyzr=%e %e %e %e, g,rho,dx=%e %e %e, slope,dP01,err=%e %e %lld",
//									i, j, k, x, y, z, r, g, rho, dx, slope, dP01, err);
						}
						dP12 = invLimiter(dP01, slope, &err);
						P2 = P1 + dP12;
						if(i == K0 && j == K0 + 2)
							TRACEF("i,j,k= %lld %lld %lld   K0,K2= %lld %lld    P0,P1,dP01= %e %e %e    slope,slope/dP01= %e %e    dP12,P2= %e %e",
									i, j, k, K0, K2, P0, P1, dP01, slope, slope / dP01, dP12, P2);
//					if(err)
//						if(i == K0 && k == K0 + 2)
//						{
//							TRACEF("ijk=%lld %lld %lld, xyzr=%e %e %e %e, g,rho,dx=%e %e %e, slope,dP01,err=%e %e %lld",
//									i, j, k, x, y, z, r, g, rho, dx, slope, dP01, err);
//						}
						if(DEBUG_HSE)
							P2 = 1e3;
						SET_HSE(i, j, k, K2, P2);

						dP01 = dP12;
						P0 = P1;
						P1 = P2;
					}
				}
			}
		TRACEF("P_boundary=%e", P2);

		if(1)
			// Integrate from the center cells out along x.
			for(size_t k = K0; k < K2; k++)
			{
				z = CELLCENTER(2, k) - SphericalGravityCenter[2];
				for(size_t j = K0; j < K2; j++)
				{
					y = CELLCENTER(1, j) - SphericalGravityCenter[1];
					P0 = totEField[ELT(K0, j, k)];
					P1 = totEField[ELT(K0 + 1, j, k)];
					dP01 = P1 - P0;
					for(size_t i = K0 + 2; i < K2; i++)
					{
						x = CELLCENTER(0, i-1) - SphericalGravityCenter[0];
						r = lenl(x, y, z);
						p->interpolateDensity(&rho, r);
						double g = -x / r * SphericalGravityGetAt(r);
						double slope = g * rho * dx;
						dP12 = invLimiter(dP01, slope, &err);
						if(err)
						{
							TRACEF("ijk=%lld %lld %lld, xyzr=%e %e %e %e, g,rho,dx=%e %e %e, slope,dP01,err=%e %e %lld",
									i, j, k, x, y, z, r, g, rho, dx, slope, dP01, err);
						}
						P2 = P1 + dP12;
						if(DEBUG_HSE)
							P2 = 1e4;
						SET_HSE(i, j, k, K2, P2);

						dP01 = dP12;
						P1 = P2;
					}
				}
			}
		TRACEF("P_boundary=%e", P2);
//Calculate the internal energy from the pressure
//		if(!DEBUG_HSE)
		for(size_t index = 0; index < FIELD_SIZE; index++)
			totEField[index] /= gammaMinusOne * rhoField[index];
#undef SET_HSE
#undef DEBUG_HSE
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

	long num, ijk[3];
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

	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField,
						NULL,
						&rhoNiField, NULL, NULL);

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
			r = lenl(xyz1, SphericalGravityCenter, GridRank);
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
