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
void printHierarchy(LevelHierarchyEntry** levelArray);
void printHierarchy0(LevelHierarchyEntry** levelArray);

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

int MHDProfileInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior)
{
	// Parameters and their defaults.
	char line[MAX_LINE_LENGTH];
	int UseMetal = FALSE;
	int useGE = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	float InternalEnergy_InitialA = 0; //[BH]
	float InternalEnergy_InitialB = 0; //[BH]

	int RefineOnStartup = FALSE;
	float RefineRegionLeft[] = { DomainLeftEdge[0], DomainLeftEdge[1], DomainLeftEdge[2] };
	float RefineRegionRight[] = { DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2] };
	int nSubgridActiveZones[] = { MetaData.TopGridDims[0], MetaData.TopGridDims[1], MetaData.TopGridDims[2] };

	char profileFileName[MAX_LINE_LENGTH];
	char profileFormat[16]; //"PLAIN" or "PAH1" or "PAH2"
	char profileType[16]; //"RADIAL"
	char radiusColumnName[32]; //"radius"
	char densityColumnName[32]; //"density"
	char temperatureColumnName[32]; //"temperature"
	float burningTemperature = 0;
	float burnedRadius = 0;
	int perturbMethod = 0;
	int perturbBeltsCount = 0;
	float perturbRaise = 1.0;
	float profileAtTime = -1;
	float dipoleMoment[3] = { 0, 0, 0 };
	float dipoleCenter[3] = { 0, 0, 0 };

	bool reinit = fptr == NULL;
	if(fptr == NULL)
		fptr = fopen("params.enzo", "r");

	*profileFileName = *profileFormat = *profileType = '\0';
	*radiusColumnName = *densityColumnName = *temperatureColumnName = '\0';

	//
	// Read initialization parameters.  Restart params are read and written
	// in ReadParameterFile.C and WriteParameterFile.
	//
	int ret = 0, ret2 = 0;
	while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
	{
		ret = 0;
//		Use PSYM or FSYM for floats, ISYM for ints.
		ret += sscanf(line, "ProfileFileName = %s", &profileFileName);
		ret += sscanf(line, "ProfileType = %s", &profileType);
		ret += sscanf(line, "ProfileFormat = %s", &profileFormat);
		ret += sscanf(line, "RadiusColumnName = %s", &radiusColumnName);
		ret += sscanf(line, "DensityColumnName = %s", &densityColumnName);
		ret += sscanf(line, "TemperatureColumnName = %s", &temperatureColumnName);
		ret += sscanf(line, "BurningTemperature = %"FSYM, &burningTemperature);
		ret += sscanf(line, "BurnedRadius = %"PSYM, &burnedRadius);
		ret += sscanf(line, "BurnedSpherePerturbationMethod = %"PSYM, &perturbMethod);
		ret += sscanf(line, "BurnedPerturbationBeltsCount = %"PSYM, &perturbBeltsCount);
		ret += sscanf(line, "BurnedPerturbationRaise = %"PSYM, &perturbRaise);
		ret += sscanf(line, "ProfileAtTime = %"FSYM, &profileAtTime);

		ret += sscanf(line, "BA = %"FSYM" %"FSYM" %"FSYM, BA, BA + 1, BA + 2);

		ret += sscanf(line, "InternalEnergy_InitialA = %"FSYM, &InternalEnergy_InitialA);	//[BH]
		ret += sscanf(line, "InternalEnergy_InitialB = %"FSYM, &InternalEnergy_InitialB);	//[BH]

		ret += sscanf(line, "MHDProfileRefineOnStartup  = %"ISYM"", &RefineOnStartup);
		ret += sscanf(line, "MHDProfileRefineAtStartupLeft  = %"PSYM" %"PSYM" %"PSYM, RefineRegionLeft,
						RefineRegionLeft + 1, RefineRegionLeft + 2);
		ret += sscanf(line, "MHDProfileRefineAtStartupRight = %"PSYM" %"PSYM" %"PSYM, RefineRegionRight,
						RefineRegionRight + 1, RefineRegionRight + 2);

		ret += sscanf(line, "MagneticDipoleMoment = %"FSYM" %"FSYM" %"FSYM, dipoleMoment, dipoleMoment + 1,
						dipoleMoment + 2);
		ret += sscanf(line, "MagneticDipoleCenter = %"FSYM" %"FSYM" %"FSYM, dipoleCenter, dipoleCenter + 1,
						dipoleCenter + 2);

//		ret2 += sscanf(line, "MHDBlastMetalDensityA = %"PSYM, &MetalDensityA);
//		ret2 += sscanf(line, "MHDBlastMetalDensityB = %"PSYM, &MetalDensityB);
//		ret2 += sscanf(line, "MHDBlastMetalOffsetInX = %"PSYM, &MetalOffsetInX);
//		if(ret2 > 0)
//		{
//			ret += ret2;
//			UseMetal = TRUE;
//			ret2 = 0;
//		}
	}

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

	bool usingVectorPotential = dipoleMoment[0] || dipoleMoment[1] || dipoleMoment[2];
	bool projectChildrenToParents = RefineOnStartup && MaximumRefinementLevel > 0;
	ret = SUCCESS;

	if(reinit){
		LevelHierarchyEntry** LevelArray = (LevelHierarchyEntry**) (Outfptr);
		LevelArrayIterator lit = LevelArrayIterator(LevelArray);
		for(grid* g = lit.firstFromTop(); g; g = lit.next())
		{
			g->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType, radiusColumnName,
										densityColumnName, temperatureColumnName, burningTemperature, burnedRadius,
										profileAtTime, dipoleMoment, dipoleCenter, usingVectorPotential);
		}

		return SUCCESS;
	}

	// Ceate a level array stub and initialize it with the TopGrid.
	LevelHierarchyEntry *levelArray[MAX_DEPTH_OF_HIERARCHY];
	arr_set(levelArray, MAX_DEPTH_OF_HIERARCHY, NULL);
	AddLevel(levelArray, &TopGrid, 0);
	LevelArrayIterator it = LevelArrayIterator(levelArray, 1);
	LevelHierarchyEntry* currentEntry = NULL;
	int currentLevel = 0;
	// The following loop initializes the TopGrid on the first pass and
	// continues if refinement is required at startup and if it is necessary.
	// After the loop, nLeves has the number of actual levels created.
	for(grid* g = it.firstFromTop(); g; g = it.next())
	{
		// A child level may not have been created during the
		// previous iteration.  If so, break out off the loop,
		// otherwise initialize the grids on this level.
		for(; g; g = it.next())
		{
			g->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType, radiusColumnName,
										densityColumnName, temperatureColumnName, burningTemperature, burnedRadius,
										profileAtTime, dipoleMoment, dipoleCenter, usingVectorPotential);
		}

		if(!RefineOnStartup || it.currentLevel >= MaximumRefinementLevel || it.currentLevel >= MAX_DEPTH_OF_HIERARCHY)
			break;

		g = it.prev();
		printHierarchy(levelArray);
		// Create a child level of parentLevel:
		if(RebuildHierarchy(&MetaData, levelArray, it.currentLevel) == FAIL)
		{
			ENZO_FAIL("Error in RebuildHierarchy.");
		}
		it.nLevels++;
		printHierarchy0(levelArray);
		CommunicationBarrier();
	}

	if(usingVectorPotential)
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
		MHD_ProjectB = !usingVectorPotential;
		MHD_ProjectE = usingVectorPotential;

		// We start this loop from the level we ended the previous loop.
		// It will not execute if numLevel==0, i.e. if the top grid hasn't been refined.
		grid* parent;
		for(grid* g = it.firstFromFinest(&parent); g; g = it.next(&parent))
		{
			if(g->ProjectSolutionToParentGrid(*parent) == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in grid->ProjectSolutionToParentGrid.");
		}
		MHD_ProjectB = origMHD_ProjectB;
		MHD_ProjectE = origMHD_ProjectE;
	}

	// Without reinement the following three loop will execute
	// once for level 0, where the TopGrid is.
	if(usingVectorPotential)
	{
		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->MHD_Curl() == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in MHD_Curl\n");

		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->CenterMagneticField() == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in CenterMagneticField\n");

		for(grid* g = it.firstFromTop(); g; g = it.next())
			if(g->MHDProfileInitializeGrid2(profileFileName, profileFormat, profileType, radiusColumnName,
											densityColumnName, temperatureColumnName, burningTemperature, burnedRadius,
											profileAtTime, dipoleMoment, dipoleCenter, usingVectorPotential) == FAIL)
				ENZO_FAIL("MHDProfileRefineOnStartup: error in MHDProfileInitializeGrid2\n");
	}

//
// Labels and Units.  (For IO.)
//
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
	for(int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
		DataUnits[i] = NULL;

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
// Done with labels and units.

//	if(UseSphericalGravity)
//	{
//		// Initialize the potential.
//		LevelHierarchyEntry tg = { NULL, TopGrid.GridData, &TopGrid };
//		LevelHierarchyEntry* topLevel[] = { &tg };
//		//ret = SphericalGravityComputePotential(topLevel, &MetaData);
//		ret = SphericalGravityComputePotential(levelArray, &MetaData);
//	}

	return ret;
}

