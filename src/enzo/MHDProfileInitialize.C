/***********************************************************************
 *
 *  Initializes the MHD problem (ProblemType 501) from external radial
 *  profiles.
 *
 *  written by: boyan hristov
 *              based on MHDBlastInititalize by David Collins
 *  date:       2018-
 *  modified1:
 *
 *  PURPOSE:  Problem initializer
 *
 *  PARAMETERS:
 *
 *	ProblemType = 501
 *       ProfileFileName = "myprofiledata.txt"
 *       ProfileFormat = "PLAIN" or "PAH1" or "PAH2"
 *       ProfileType = "RADIAL"
 *       RadiusColumnName = "radius"
 *       DensityColumnName = "density"
 *       TemperatureColumnName = "temperature"
 *       BurnedRadius = 1e5
 *       BurningTemperature = 8.5e9
 *
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelArrayIterator.h"
#include "DebugTools.h"

class ExternalBoundary;

/* This initializes genearlized discontinuities.  Good for blast waves, Kelving Helmholtz, shock tubes.  Probably others.

 Code flow:
 1.) Declare/Define Defaults for  parameters.
 2.) Read parameters from file.
 3.) Calculate TotalEnergy from quantity given (options are Total Energy, Gas Energy, Pressure.)
 4.) Set up data labels, units.
 5.) Declare Hierarchy Object.
 6.) Define linked list
 7.) Call initializer on each level.
 8.) Project to parent.

 */

//in MHD_ObliqueRoutines.C
//void RotateVector(float * Vector, float * Normal);
//int SetupNormal(float Normal[], float MHDBlastCenter[3], TopGridData & MetaData);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int CommunicationBarrier();
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level);
int MHDCTSetupFieldLabels();
int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData);

int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z)
{
	*Bx = BA[0];
	*By = BA[1];
	*Bz = BA[2];
	return SUCCESS;
}

int MHDProfileInitExactB(float B[3], FLOAT x[3])
{
	return MHDProfileInitExactB(B, B + 1, B + 2, x[0], x[1], x[2]);
}

int MHDProfileSetLabels(bool useGE, bool UseMetal)
{
	char *DensName = "Density";
	char *TEName = "TotalEnergy";
	char *Vel1Name = "x-velocity";
	char *Vel2Name = "y-velocity";
	char *Vel3Name = "z-velocity";
	char *GPotName = "GravPotential";
	char *BxName = "Bx";
	char *ByName = "By";
	char *BzName = "Bz";
	char *PhiName = "Phi";
	char *GEName = "GasEnergy";
	char *MetalName = "Metal_Density";
	char *Density_56NiName = "Density_56Ni";   //[BH]
	char *MHDCT_BxName = "BxF";
	char *MHDCT_ByName = "ByF";
	char *MHDCT_BzName = "BzF";
	char *MHDCT_ExName = "Ex";
	char *MHDCT_EyName = "Ey";
	char *MHDCT_EzName = "Ez";

	// Ignore the unit strings.
	arr_set(DataUnits, MAX_NUMBER_OF_BARYON_FIELDS, NULL);

	int i = 0;
	DataLabel[i++] = DensName; //"g/cm**3";

	if(EquationOfState == 0)
		DataLabel[i++] = TEName; //"erg/g";

	if(useGE)
		DataLabel[i++] = GEName; //"erg/g";

	DataLabel[i++] = Vel1Name; //"cm/s";
	DataLabel[i++] = Vel2Name; //"cm/s";
	DataLabel[i++] = Vel3Name; //"cm/s";

	if(UseMetal)
		DataLabel[i++] = MetalName;

	if(WritePotential)
		DataLabel[i++] = GPotName;

	if(UseBurning)
		DataLabel[i++] = Density_56NiName; //[BH]

	if(UseMHD)
	{
		DataLabel[i++] = BxName;
		DataLabel[i++] = ByName;
		DataLabel[i++] = BzName;
		if(HydroMethod == MHD_RK)
			DataLabel[i++] = PhiName;
	}
	if(UseMHDCT)
		MHDCTSetupFieldLabels();

	return SUCCESS;
}

int MHDProfileInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior, bool BeforeGridDistribution)
{
	// Parameters and their defaults.
	char line[MAX_LINE_LENGTH];
	int UseMetal = FALSE;
	int useGE = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	char ProfileFileName[MAX_LINE_LENGTH];
	char ProfileFormat[16]; //"PLAIN" or "PAH1" or "PAH2"
	char ProfileType[16]; //"RADIAL"
	char RadiusColumnName[32]; //"radius"
	char DensityColumnName[32]; //"density"
	char InternalEnergyColumnName[32];
	char TemperatureColumnName[32]; //"temperature"
	float BurningTemperature = 0;
	int PerturbMethod = 0;
	float ProfileAtTime = -1;
	float dipoleMoment[3] = { 0, 0, 0 };
	float dipoleCenter[3] = { 0, 0, 0 };

//	bool reinit = fptr == NULL;
//	if(fptr == NULL)
//		fptr = fopen("params.enzo", "r");

	*ProfileFileName = *ProfileFormat = *ProfileType = '\0';
	*RadiusColumnName = *DensityColumnName = *TemperatureColumnName = '\0';
	*InternalEnergyColumnName = '\0';

	int ret = 0, ret2 = 0;
	while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
	{
		ret = 0;
//		Use PSYM or FSYM for floats, ISYM for ints.
		ret += sscanf(line, "ProfileFileName = %s", &ProfileFileName);
		ret += sscanf(line, "ProfileType = %s", &ProfileType);
		ret += sscanf(line, "ProfileFormat = %s", &ProfileFormat);
		ret += sscanf(line, "RadiusColumnName = %s", &RadiusColumnName);
		ret += sscanf(line, "DensityColumnName = %s", &DensityColumnName);
		ret += sscanf(line, "InternalEnergyColumnName = %s", &InternalEnergyColumnName);
		ret += sscanf(line, "TemperatureColumnName = %s", &TemperatureColumnName);
		ret += sscanf(line, "BurningTemperature = %"FSYM, &BurningTemperature);
		ret += sscanf(line, "ProfileAtTime = %"FSYM, &ProfileAtTime);

		ret += sscanf(line, "BA = %"FSYM" %"FSYM" %"FSYM, BA, BA + 1, BA + 2);

		ret += sscanf(line, "MagneticDipoleMoment = %"FSYM" %"FSYM" %"FSYM, dipoleMoment, dipoleMoment + 1,
						dipoleMoment + 2);
		ret += sscanf(line, "MagneticDipoleCenter = %"FSYM" %"FSYM" %"FSYM, dipoleCenter, dipoleCenter + 1,
						dipoleCenter + 2);
	}

	MHDProfileSetLabels(useGE, UseMetal);
	if(ParallelRootGridIO && BeforeGridDistribution)
	{
		return TopGrid.GridData->MHDProfileInitializeGrid(NULL, BurningTemperature, InitialBurnedRadius, dipoleMoment,
															dipoleCenter, InitBWithVectorPotential);
	}

	MHDInitialProfile p = MHDInitialProfile();
	p.init(RadiusColumnName, DensityColumnName, InternalEnergyColumnName);
	p.read(ProfileFileName, ProfileFormat, ProfileAtTime);
	printf("Profile: %lld data rows in %lld columns.\n", p.nRows, p.nCols);

	if(BurningDiffusionRateReduced <= 0)
		BurningDiffusionRateReduced = BurningDiffusionRate / (DomainRightEdge[0] - DomainLeftEdge[0])
				* MetaData.TopGridDims[0];
	if(BurningReactionRateReduced <= 0)
		BurningReactionRateReduced = BurningReactionRate * (DomainRightEdge[0] - DomainLeftEdge[0])
				/ MetaData.TopGridDims[0];

	if(SphericalGravityOuterRadius <= 0)
	{
		// If outer radius is not specified, set it to the minimal
		// radius to cover the entire domain, i.e. the distance
		// from the gravity center to the farthest domain vertice.
		SphericalGravityOuterRadius = distMaxl(DomainLeftEdge, DomainRightEdge, SphericalGravityCenter,
												MetaData.TopGridRank);
	}

	InitBWithVectorPotential = InitBWithVectorPotential && (dipoleMoment[0] || dipoleMoment[1] || dipoleMoment[2]);
	bool projectChildrenToParents = RefineOnStartup && MaximumRefinementLevel > 0;
	int maxRefLevel = min(MaximumRefinementLevel, MAX_DEPTH_OF_HIERARCHY);
	maxRefLevel = (maxRefLevel >= 0) ? maxRefLevel : 2;
	maxRefLevel *= (bool) RefineOnStartup;

	ret = SUCCESS;

	RebuildHierarchyIterator rhit = RebuildHierarchyIterator(maxRefLevel, &TopGrid, &MetaData);
	for(grid* g = rhit.first(); g; g = rhit.next())
	{
		if(rhit.startingNewLevel && rhit.currentLevel == 0)
			PRINT_HIERARCHY0(rhit.levelArray);

		g->MHDProfileInitializeGrid(&p, BurningTemperature, InitialBurnedRadius, dipoleMoment, dipoleCenter,
									InitBWithVectorPotential);
	}

	TRACEF("HIERARCHY BUILT");

	LevelArrayIterator it = LevelArrayIterator(rhit.levelArray);
	if(UseSphericalGravity)
	{
		// Initialize the potential.
		ret = SphericalGravityComputePotential(rhit.levelArray, &MetaData);
	}

	if(InitBWithVectorPotential)
	{
		for(grid* g = it.firstFromTop(); g; g = it.next())
			g->InitializeMagneticUniformFieldVectorPotential(BA, 0);

		for(grid* g = it.firstFromTop(); g; g = it.next())
			g->InitializeMagneticDipoleVectorPotential(dipoleMoment, dipoleCenter, 1);
	}

	if(projectChildrenToParents)
	{
		// Restore the consistency among levels by projecting each layer to its parent,
		// starting from the finest level.
		// We want to project the vector potential and take the curl afterwards.
		// Store the original project flags for the MHD fields and set them so
		// that the vector potential gets projected.
		int origMHD_ProjectB = MHD_ProjectB;
		int origMHD_ProjectE = MHD_ProjectE;
		MHD_ProjectB = !InitBWithVectorPotential;
		MHD_ProjectE = InitBWithVectorPotential;

		// We start this loop from the level we ended the previous loop.
		// It will not execute if numLevel==0, i.e. if the top grid hasn't been refined.
		grid* parent;
		for(grid* g = it.firstFromFinest(&parent); g; g = it.next(&parent))
		{
			if(parent)
				if(g->ProjectSolutionToParentGrid(*parent) == FAIL)
					ENZO_FAIL("MHDProfileRefineOnStartup: error in grid->ProjectSolutionToParentGrid.");
		}
		MHD_ProjectB = origMHD_ProjectB;
		MHD_ProjectE = origMHD_ProjectE;
	}

// Without reinement the following three loop will execute
// once for level 0, where the TopGrid is.
	if(InitBWithVectorPotential)
	{
		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->MHD_Curl() == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in MHD_Curl\n");

		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->CenterMagneticField() == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in CenterMagneticField\n");
	}

	if(InitBWithVectorPotential
			|| (UseSphericalGravity && InitRadialPressureFromCentral && *InternalEnergyColumnName == '\0'))
	{
		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->MHDProfileInitializeGrid2(&p, BurningTemperature, InitialBurnedRadius, dipoleMoment, dipoleCenter,
											InitBWithVectorPotential) == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in MHDProfileInitializeGrid2\n");
	}

	p.free();
	return ret;
}

void WriteInitialProfile(char* name, FLOAT* RR, FLOAT* RHO, FLOAT* GG, FLOAT* PP, FLOAT* UU, size_t n, FLOAT K,
FLOAT gamma)
{
	if(MyProcessorNumber != ROOT_PROCESSOR)
		return;

	const int M = 10;
	size_t num;
	FLOAT r;
	float rho, P, U, g;
	char filename[FILENAME_MAX];
	sprintf(filename, "%s.initialProfile", name);
	FILE* file = fopen(filename, "w");
	if(file == NULL)
	{
		ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
	}
	fprintf(file, "# rownum   r   rho   P   U   g\n");
	for(size_t i = 0; i < n; i++)
	{
		num = i + 1;
		r = RR[i];
		rho = RHO[i];
		P = PP[i];
		U = UU[i];
		g = GG[i];
		fprintf(file, "%lld   %12.9e  %12.9e  %12.9e  %12.9e  %12.9e\n", num, r, rho, P, U, g);
	}
	if(file != stderr && file != stdout)
	{
		fclose(file);
		fprintf(stderr, "'%s' written.\n", filename);
	}

	sprintf(filename, "%s.initialProfile.py", name);
	file = fopen(filename, "w");
	if(file == NULL)
	{
		ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
	}

	fprintf(file, "type('', (), dict(\n");
	fprintf(file, "time = %f,\n", 0);
	fprintf(file, "gamma = %f,\n", gamma);
	fprintf(file, "K = %12.9e,\n", K);
	fprintf(file, "cols = 'r rho P U g'.split(),\n");

	fprintf(file, "%s = np.array([", "r");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", RR[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "rho");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", RHO[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "P");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", PP[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "U");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", UU[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "g");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", GG[i]);
	}
	fprintf(file, "\n])\n))\n");

	if(file != stderr && file != stdout)
	{
		fclose(file);
		fprintf(stderr, "'%s' written.\n", filename);
	}

}

//int SphericalGravityIntegratePressure(float* pressure_out, float* gasE_out, FLOAT* radius, float* gravity,
//	float* density, size_t n,
//	float startingPressure, bool startAtCenter)
//{
//	float P_c = 0;
//	float rho_c = 0;
//	FLOAT xyz[MAX_DIMENSION];
//	arr_cpy(xyz, SphericalGravityCenter, MAX_DIMENSION);
//
//	if(InitRadialPressureFromCentral > 0)
//	{
//		P_c = InitRadialPressureFromCentral;
//		if(MyProcessorNumber == ROOT_PROCESSOR)
//			printf("P_c = %e (parameter)", P_c);
//	}
//	else
//	{
////			size_t i_c = getCellIndex(xyz);
////			float rho_c = rhoField[i_c];
//		FLOAT rho_c;
//		retcode = p->interpolateDensity(&rho_c, 0); //g/cm**3
//		P_c = EOSPolytropicFactor * POW(rho_c, Gamma);
//		if(MyProcessorNumber == ROOT_PROCESSOR)
//			printf("P_c = %e (polytropic, = K*rho**gamma, K=%e, rho=%e, gamma=%5.3f)", P_c, EOSPolytropicFactor,
//					rho_c, gamma);
//	}
//
//	FLOAT* P_r = new float[SphericalGravityActualNumberOfBins];
//	FLOAT* U_r = new float[SphericalGravityActualNumberOfBins];
//	FLOAT* rr = NULL;
//	FLOAT* dd = NULL;
//	FLOAT* gg = NULL;
//	if(ProcessorNumber == ROOT_PROCESSOR && ID == 0 && CellWidth[0][0] == TopGridDx[0])
//	{
//		rr = new FLOAT[SphericalGravityActualNumberOfBins];
//		gg = new FLOAT[SphericalGravityActualNumberOfBins];
//		dd = new FLOAT[SphericalGravityActualNumberOfBins];
//	}
//
//	FLOAT P2, P1, P, r2, r1, r, U2, U1, U; // r2 < r1 < r
//	FLOAT b, a, bb, aa, g, rho, g_rho; // b:=r-r2 > a:=r-r1
//
//	r2 = 0; // Start from the center.
//	P_r[0] = P2 = 0; // Offfset the pressure later.
//
//	//xyz[0] = r2 + SphericalGravityCenter[0];
//	//size_t index = getCellIndex(xyz);
//	retcode = p->interpolateDensity(&rho, r2); //g/cm**3
//	g_rho = g = SphericalGravityGetAt(r2);
////		TRACEF("r[0], P[0], g[0]rho[0], g[0], rho[0] = %e, %e, %e = %e * %e", r2, P2, g_rho, g, rho);
//	if(rr)
//	{
//		rr[0] = r2;
//		gg[0] = g;
//		dd[0] = rho;
//	}
//
//	r1 = SphericalGravityBinLeftEdges[1];
////		xyz[0] = r1 + SphericalGravityCenter[0];
////		index = getCellIndex(xyz);
////		rho1 = rhoField[index];
//	retcode = p->interpolateDensity(&rho, r1); //g/cm**3
//	g = SphericalGravityGetAt(r1);
//	g_rho = g * rho;
//	if(rr)
//	{
//		rr[1] = r1;
//		gg[1] = g;
//		dd[1] = rho;
//	}
//
//	a = r1 - r;
//	P1 = P2 - 0.5 * a * g * g_rho;
//	P1 -= 2. / 3. * M_PI * SphericalGravityConstant * a * a * rho * rho;
//
//	P_r[1] = P1;
////		TRACEF("r[1], P[1], dP[1], g, rho = %e, %e, %e = %e * %e", r1, P1, P1 - P2, g, rho);
//
//	for(size_t i = 2; i < SphericalGravityActualNumberOfBins; i++)
//	{
//		r = SphericalGravityBinLeftEdges[i];
//		b = r - r2;
//		a = r - r1;
//		bb = square(b);
//		aa = square(a);
//
////			xyz[0] = r0 + SphericalGravityCenter[0];
////			index = getCellIndex(xyz);
////			rho = rhoField[index];
//		retcode = p->interpolateDensity(&rho, r); //g/cm**3
//		g = -SphericalGravityGetAt(r);
//		g_rho = g * rho;
//		P = (bb * (P1 + a * g_rho) - aa * (P2 + b * g_rho)) / (bb - aa);
////			TRACEF("%lld: P2 P1 P = %e %e %e ; r2 r1 r = %e %e %e ; g rho = %e %e", i, P2, P1, P, r2, r1, r, g, rho);
//
//		if(rr)
//		{
//			rr[i] = r;
//			gg[i] = g;
//			dd[i] = rho;
//		}
//
//		P_r[i] = P;
//		r2 = r1;
//		r1 = r;
//		P2 = P1;
//		P1 = P;
//	}
//
//	// P and rho have the outermost values from the P integration loop.
//	// Caclculate correction as if intergrating from the surface inwards.
//	const double K53 = 3.161128e+12;
//	const double K43 = 4.881986e+14;
//	double K = 3.828648955e+14;
//	P1 = K * POW(rho, Gamma);
//	P = P1 - P2; // add to P_r to get the right boundary conditions.
//	P = K * POW(2e9, Gamma);
////		TRACEF("Pressure shift: %e = %e - %e", P2, P1, P);
//	// Shift pressure and replace it with the gas energy.
//	for(size_t i = 0; i < SphericalGravityActualNumberOfBins; i++)
//	{
//		a = P_r[i] + P;
//		r = SphericalGravityBinLeftEdges[i];
//		retcode = p->interpolateDensity(&rho, r); //g/cm**3
//		b = a / rho / gammaMinusOne;
//		P_r[i] = a;
//		U_r[i] = b;
//	}
//
//	if(rr)
//	{
//		WriteInitialProfile("initialProfile", rr, dd, gg, P_r, U_r, SphericalGravityActualNumberOfBins, K,
//							(double)Gamma);
//	}
//	for(int k = 0; k < GridDimension[2]; k++)
//	{
//		FLOAT z = CELLCENTER(2, k);
//		FLOAT zz = square(z - SphericalGravityCenter[2]);
//		for(int j = 0; j < GridDimension[1]; j++)
//		{
//			FLOAT y = CELLCENTER(1, j);
//			FLOAT yy_zz = square(y - SphericalGravityCenter[1]) + zz;
//			for(int i = 0; i <= GridDimension[0]; i++)
//			{
//				FLOAT x = CELLCENTER(0, i);
//				FLOAT r = sqrt(square(x - SphericalGravityCenter[0]) + yy_zz);
//				int index = ELT(i, j, k);
//
//				float rho = rhoField[index];
//				int retval;
//				if(p->internalEnergyData)
//				{
//					retval = p->interpolateInternalEnergy(&U, r);
//				}
//				else
//				{
//					U = 0;
////						size_t rbin = SphericalGravityComputeBinIndex(r);
////						r1 = SphericalGravityBinLeftEdges[rbin];
////						r2 = SphericalGravityBinLeftEdges[rbin + 1];
////						// P_r[] has the specific gas energy now.
////						U1 = U_r[rbin];
////						U2 = U_r[rbin + 1];
////						U = ((r2 - r) * U1 + (r - r1) * U2) / (r2 - r1);
//				}
//				float gasE = U; // / rho / gammaMinusOne;
//				totEField[index] = U; // add kinetic and magnetic energy below
//			}
//		}
//	}
//
//	if(debug)
//	{
//		for(size_t i = 0; i < SphericalGravityActualNumberOfBins; i++)
//		{
//			r = SphericalGravityBinLeftEdges[i];
//			retcode = p->interpolateDensity(&rho, r); //g/cm**3
//			P = P_r[i];
//			U = U_r[i];
//			TRACEGF0("Intergrated ressure: %lld: r=%e P=%e U=%e rho=%e gamma-1=%f", i, r, P, U, rho, gammaMinusOne);
//		}
//	}
//
//	delete P_r;
//	delete U_r;
//	delete rr;
//	delete gg;
//	delete dd;
//
//}
//
