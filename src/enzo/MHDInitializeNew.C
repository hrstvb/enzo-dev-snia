/***********************************************************************
 /
 /  INITIALIZE A NEW SIMULATION
 /
 /  written by: Greg Bryan
 /  date:       November, 1994
 /  modified1:  Robert Harkness
 /  date:       September 2004
 /  modified2:  Stephen Skory
 /  date:       May, 2008
 /  modified3:  Alexei Kritsuk
 /  date:       May, 2008
 /
 /  PURPOSE:
 /
 /  RETURNS: SUCCESS or FAIL
 /
 ************************************************************************/

// This routine intializes a new simulation based on the parameter file.
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>

#include "DebugMacros.h"
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
#include "CommunicationUtilities.h"

// Function prototypes

int InitializeMovieFile(TopGridData &MetaData, HierarchyEntry &TopGrid);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData, char *Filename = NULL);
void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid);
int SetDefaultGlobalValues(TopGridData &MetaData);
//int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);
int CommunicationBroadcastValue(PINT *Value, int BroadcastProcessor);
int TracerParticleCreation(FILE *fptr, HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHDCT_ParameterJuggle(); //updates old style MHDCT parameter files to reflect new values
void PrintMemoryUsage(char *str);

int InitializeNew_pass1(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr);
int InitializeNew_pass2(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr);
int InitializeNew_TopGrid(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr);
int InitializeNew_Exterior(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr);
int InitializeNew_PartitionGrids(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior, float *Initialdt, FILE *fptr, FILE *Outfptr, bool partitionNestedGrids);

// Character strings

char outfilename[] = "amr.out";

int InitializeNew(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt)
{
	TRACEF("INITIALIZING NEW PROBLEM  BEGIN ...");
	FILE *fptr, *Outfptr;
	int dim, i;

	// Open parameter file
	if((fptr = fopen(filename, "r")) == NULL)
		ENZO_VFAIL("Error opening parameter file '%s'.", filename);

	// Clear OutputLog
	FILE *sptr;
	if(MyProcessorNumber == ROOT_PROCESSOR)
	{
		sptr = fopen("OutputLog", "w");
		fclose(sptr);
	}

	// Open output file
	if(MyProcessorNumber == ROOT_PROCESSOR)
		if((Outfptr = fopen(outfilename, "w")) == NULL)
			ENZO_VFAIL("Error opening parameter output file %s\n", outfilename)

	// set the default MetaData values
	SetDefaultGlobalValues(MetaData);

	// Read the MetaData/global values from the Parameter file
	if(ReadParameterFile(fptr, MetaData, Initialdt) == FAIL)
		ENZO_FAIL("Error in ReadParameterFile.");

	//Ensure old style MHD_CT parameter files still work.
	if(MHDCT_ParameterJuggle() == FAIL)
		ENZO_FAIL("Invalid parameter from old style MHD CT");

	// Set the number of particle attributes, if left unset
	if(NumberOfParticleAttributes == INT_UNDEFINED || NumberOfParticleAttributes == 0)
	{
		if(StarParticleCreation || StarParticleFeedback)
		{
			NumberOfParticleAttributes = 3;
			if(StarMakerTypeIaSNe)
				NumberOfParticleAttributes++;
			if(StarMakerTypeIISNeMetalField)
				NumberOfParticleAttributes++;
		}
		else
		{
			NumberOfParticleAttributes = 0;
		}
	}

	// Give unset parameters their default values
	for(dim = 0; dim < MAX_DIMENSION; dim++)
	{
		if(RefineRegionLeftEdge[dim] == FLOAT_UNDEFINED)
			RefineRegionLeftEdge[dim] = DomainLeftEdge[dim];
		if(RefineRegionRightEdge[dim] == FLOAT_UNDEFINED)
			RefineRegionRightEdge[dim] = DomainRightEdge[dim];
	}

	// If the problem reads in a restart dump, then skip over the following
	if(ProblemType != 40 && ProblemType != 51)
		InitializeNew_TopGrid(filename, TopGrid, MetaData, Exterior, Initialdt, fptr, Outfptr);

	// Call a specific problem initializer
	PrintMemoryUsage("Call problem init");
	int ret = InitializeNew_pass1(filename, TopGrid, MetaData, Exterior, Initialdt, fptr, Outfptr);

	if(ret == INT_UNDEFINED)
		ENZO_VFAIL("Problem Type %"ISYM" undefined.\n", ProblemType)

	if(ret == FAIL)
		ENZO_FAIL("Error in problem initialization.");

	/* Do some error checking */
	if(MetaData.StopTime == FLOAT_UNDEFINED && MetaData.StopCycle == INT_UNDEFINED)
		ENZO_FAIL("StopTime nor StopCycle ever set.");
	if(MetaData.StopCycle != INT_UNDEFINED && MetaData.StopTime == FLOAT_UNDEFINED)
		MetaData.StopTime = huge_number;

	int nFields = TopGrid.GridData->ReturnNumberOfBaryonFields();
	if(nFields >= MAX_NUMBER_OF_BARYON_FIELDS)
	{
		ENZO_VFAIL("NumberOfBaryonFields (%"ISYM") + 1 exceeds " "MAX_NUMBER_OF_BARYON_FIELDS (%"ISYM").\n", nFields,
					MAX_NUMBER_OF_BARYON_FIELDS)
	}

	PrintMemoryUsage("1st Initialization done");
	InitializeNew_Exterior(filename, TopGrid, MetaData, Exterior, Initialdt, fptr, Outfptr);
	PrintMemoryUsage("Exterior set");

	// Set values that were left undefined (above)
	if(MetaData.TimeLastDataDump == FLOAT_UNDEFINED)
		MetaData.TimeLastDataDump = MetaData.Time - MetaData.dtDataDump * 1.00001;
	if(MetaData.TimeLastHistoryDump == FLOAT_UNDEFINED)
		MetaData.TimeLastHistoryDump = MetaData.Time - MetaData.dtHistoryDump;
	if(MetaData.TimeLastTracerParticleDump == FLOAT_UNDEFINED)
		MetaData.TimeLastTracerParticleDump = MetaData.Time - MetaData.dtTracerParticleDump;
	if(MetaData.CycleLastDataDump == INT_UNDEFINED)
		MetaData.CycleLastDataDump = MetaData.CycleNumber - MetaData.CycleSkipDataDump;
	if(MetaData.CycleLastHistoryDump == INT_UNDEFINED)
		MetaData.CycleLastHistoryDump = MetaData.CycleNumber - MetaData.CycleSkipHistoryDump;

	// Make changes required for Zeus solver, and turn the TotalEnergy
	// variable (should be renamed just Energy) into GasEnergy
	if(HydroMethod == Zeus_Hydro && ProblemType != 10 &&  // BWO (Rotating cylinder)
			ProblemType != 11 &&  // BWO (radiating shock)
			ProblemType != 13 &&  // BWO (Rotating Sphere)
			ProblemType != 20 && ProblemType != 27 && ProblemType != 30 && ProblemType != 31 && // BWO (isolated galaxies)
			ProblemType != 60 && ProblemType != 106 && //AK
			ProblemType != 108 && //Yuan (Cluster)
			ProblemType != 501 //BH (MHD WD)
					)
		ConvertTotalEnergyToGasEnergy(&TopGrid);

	// If using StarParticles, set the number to zero
	// (assuming it hasn't already been set)
	if(NumberOfStarParticles != 0)
		if(StarParticleCreation || StarParticleFeedback)
			NumberOfStarParticles = 0;

	// Convert minimum initial overdensity for refinement to mass
	// (unless MinimumMass itself was actually set)
	for(i = 0; i < MAX_FLAGGING_METHODS; i++)
		if(MinimumMassForRefinement[i] == FLOAT_UNDEFINED)
		{
			MinimumMassForRefinement[i] = MinimumOverDensityForRefinement[i];
			for(dim = 0; dim < MetaData.TopGridRank; dim++)
				MinimumMassForRefinement[i] *= (DomainRightEdge[dim] - DomainLeftEdge[dim]) /
				float(MetaData.TopGridDims[dim]);
		}

	// Check for the creation of tracer particles
	// Tracer particles will not be created at this point if ||rgio in ON
	if(TracerParticleCreation(fptr, TopGrid, MetaData) == FAIL)
		ENZO_FAIL("Error in TracerParticleCreation");

	// Write the MetaData/global values to the Parameter file
	if(MyProcessorNumber == ROOT_PROCESSOR)
		if(WriteParameterFile(Outfptr, MetaData) == FAIL)
		{
			ENZO_FAIL("Error in WriteParameterFile.");
		}

	if(debug)
		printf("InitializeNew: Initial grid hierarchy set\n");

//	// Walk the grids
//	HierarchyEntry *CurrentGrid;
//	int gridcounter = 0;
//
//	CurrentGrid = &TopGrid;
//	while(CurrentGrid != NULL)
//	{
////		if(debug)
//		printf("InitializeNew: Partition Initial Grid %"ISYM"\n", gridcounter);
//
//		if(CurrentGrid->NextGridThisLevel == NULL)
//			CommunicationPartitionGrid(CurrentGrid, gridcounter);
//
//		gridcounter++;
//
//		if(PartitionNestedGrids)
//			CurrentGrid = CurrentGrid->NextGridNextLevel;
//		else
//			CurrentGrid = NULL;
//	}
	InitializeNew_PartitionGrids(filename, TopGrid, MetaData, Exterior, Initialdt, fptr, Outfptr, (bool)PartitionNestedGrids);

	PrintMemoryUsage("Before 2nd pass");
	InitializeNew_pass2(filename, TopGrid, MetaData, Exterior, Initialdt, fptr, Outfptr);

	// Close parameter files
	fclose(fptr);
	if(MyProcessorNumber == ROOT_PROCESSOR)
		fclose(Outfptr);

	MetaData.FirstTimestepAfterRestart = FALSE;

	PrintMemoryUsage("Exit X_Init");

	CommunicationBarrier();

	// 2006-12-11 Skory bug fix for star particle miscounts
	// Added the following line:
	CommunicationBroadcastValue(&MetaData.NumberOfParticles, ROOT_PROCESSOR);

	/* If requested, initialize streaming data files. */
	InitializeMovieFile(MetaData, TopGrid);

	if(debug)
		printf("InitializeNew: Finished problem initialization.\n");

	return SUCCESS;
}

void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid)
{
	if(Grid != NULL)
	{

		Grid->GridData->ConvertTotalEnergyToGasEnergy();
		ConvertTotalEnergyToGasEnergy(Grid->NextGridThisLevel);
		ConvertTotalEnergyToGasEnergy(Grid->NextGridNextLevel);
	}
}
