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
int MHDCTSetupFieldLabels();

int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z)
{
	*Bx = *By = *Bz = 0;
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
	float burningRadius = 0;
	float profileAtTime = -1;

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
		ret += sscanf(line, "BurningRadius = %"PSYM, &burningRadius);
		ret += sscanf(line, "ProfileAtTime = %"FSYM, &profileAtTime);

		ret += sscanf(line, "InternalEnergy_InitialA = %"FSYM, &InternalEnergy_InitialA);	//[BH]
		ret += sscanf(line, "InternalEnergy_InitialB = %"FSYM, &InternalEnergy_InitialB);	//[BH]

		ret += sscanf(line, "MHDProfileRefineOnStartup  = %"ISYM"", &RefineOnStartup);
		ret += sscanf(line, "MHDProfileRefineAtStartupLeft  = %"PSYM" %"PSYM" %"PSYM, RefineRegionLeft,
						RefineRegionLeft + 1, RefineRegionLeft + 2);
		ret += sscanf(line, "MHDProfileRefineAtStartupRight = %"PSYM" %"PSYM" %"PSYM, RefineRegionRight,
						RefineRegionRight + 1, RefineRegionRight + 2);

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

	// Initialize the top grid.
	if(TopGrid.GridData->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType, radiusColumnName,
													densityColumnName, temperatureColumnName, burningTemperature,
													burningRadius, profileAtTime) == FAIL)
		ENZO_FAIL("MHDProfileInitialize: Error in MHDProfileInitializeGrid.");

	if(RefineOnStartup && MaximumRefinementLevel > 0)
	{
		// Create a new hierarchy stub with single entry at each level.
		HierarchyEntry ** Subgrid = new HierarchyEntry*[MaximumRefinementLevel];
		for(int lev = 0; lev < MaximumRefinementLevel; lev++)
			Subgrid[lev] = new HierarchyEntry;

		for(int lev = 0; lev < MaximumRefinementLevel; lev++)
			Subgrid[lev]->NextGridThisLevel = NULL;

		TopGrid.NextGridNextLevel = Subgrid[0];
		Subgrid[0]->ParentGrid = &TopGrid;
		for(int lev = 1; lev < MaximumRefinementLevel; lev++)
		{
			Subgrid[lev - 1]->NextGridNextLevel = Subgrid[lev];
			Subgrid[lev]->ParentGrid = Subgrid[lev - 1];
		}
		Subgrid[MaximumRefinementLevel - 1]->NextGridNextLevel = NULL;

		// Align the refinement region with the top grid cells and
		// calculate the number of enclosed (active) zones.
		int SubgridDims[3];
		for(int dim = 0; dim < 3; dim++)
		{
			const FLOAT& L = DomainLeftEdge[dim];
			const FLOAT& R = DomainRightEdge[dim];
			const FLOAT DX = (R - L) / MetaData.TopGridDims[dim];
			int n;

			//nint = nearest int
			n = nint((RefineRegionLeft[dim] - L) / DX);
			RefineRegionLeft[dim] = max(n * DX + L, 0);
			n = nint((RefineRegionRight[dim] - L) / DX);
			RefineRegionRight[dim] = min(n * DX + L, R);
			nSubgridActiveZones[dim] = nint((RefineRegionRight[dim] - RefineRegionLeft[dim]) / DX);
		}

		fprintf(stderr, "Subgrid Left %"GSYM" %"GSYM" %"GSYM"\n", RefineRegionLeft[0], RefineRegionLeft[1],
				RefineRegionLeft[2]);
		fprintf(stderr, "Subgrid Right %"GSYM" %"GSYM" %"GSYM"\n", RefineRegionRight[0], RefineRegionRight[1],
				RefineRegionRight[2]);
		fprintf(stderr, "nCells %"ISYM" %"ISYM" %"ISYM"\n", nSubgridActiveZones[0], nSubgridActiveZones[1],
				nSubgridActiveZones[2]);

		// Instantiate one child per parent
		for(int lev = 0; lev < MaximumRefinementLevel; lev++)
		{
			//Calculate number of cells on this level.
			for(int dim = 0; dim < MetaData.TopGridRank; dim++)
				nSubgridActiveZones[dim] *= RefineBy;

			fprintf(stderr, "MHDProfileInitialize: Refinement level=%"ISYM
			", NumberOfSubgridZones[0] = %"ISYM"\n",
					lev + 1, nSubgridActiveZones[0]);

			if(nSubgridActiveZones[0] <= 0)
			{
				printf("SedovBlast: single grid start-up.\n");
				continue;
			}

			//  Compute the full dimensions = active + ghost
			for(int dim = 0; dim < MetaData.TopGridRank; dim++)
				SubgridDims[dim] = nSubgridActiveZones[dim] + 2 * NumberOfGhostZones;

			// Create a new subgrid and initialize it
			Subgrid[lev]->GridData = new grid;
			Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
			Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims, RefineRegionLeft, RefineRegionRight,
												0);

			if(Subgrid[lev]->GridData->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType,
																radiusColumnName, densityColumnName,
																temperatureColumnName, burningTemperature,
																burningRadius, profileAtTime) == FAIL)
				ENZO_FAIL("MHDBlastInitialize: Error in MHDBlastInitializeGrid.");
		} // level

		// Make sure each grid has the best data with respect to the finer grids.
		// This projection juggle is to ensure that, regardless of how the hierarchy is evolved,
		// the field gets projected properly here.
		int MHD_ProjectEtmp = MHD_ProjectE;
		int MHD_ProjectBtmp = MHD_ProjectB;
		MHD_ProjectE = FALSE;
		MHD_ProjectB = TRUE;

		for(int lev = MaximumRefinementLevel - 1; lev > 0; lev--)
			if(Subgrid[lev]->GridData->ProjectSolutionToParentGrid(*(Subgrid[lev - 1]->GridData)) == FAIL)
				ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
		// Project to the root grid
		if(Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData)) == FAIL)
			ENZO_FAIL("Error in ProjectSolutionToParentGrid.");

		// Restore projection options to their initial.
		MHD_ProjectE = MHD_ProjectEtmp;
		MHD_ProjectB = MHD_ProjectBtmp;

	} //RefineOnStartup

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
		DataLabel[i++] = Density_56NiName;    //[BH]

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

