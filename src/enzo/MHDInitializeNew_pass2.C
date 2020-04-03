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

// Initialization function prototypes

#ifdef TRANSFER
int PhotonTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData, bool Reinitialize =
		false);
#endif /* TRANSFER */

int TurbulenceInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
int SetBaryonFields);

int DrivenFlowInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
int SetBaryonFields);

int MHD2DTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData, int SetBaryonFields);
int CollapseMHD3DInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
int SetBaryonFields);
int MHDTurbulenceInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
int SetBaryonFields);
int MHDDecayingRandomFieldInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
int SetBaryonFields);

int MHDProfileInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior, bool BeforeGridDistribution);

int CosmologySimulationReInitialize(HierarchyEntry *TopGrid, TopGridData &MetaData);
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid, TopGridData &MetaData);
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid, TopGridData &MetaData);
void PrintMemoryUsage(char *str);

int InitializeNew_pass2(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr)
{
	TRACE;
	// For problem 30, using ParallelGridIO,
	// read in data only after partitioning the grid
	if(ParallelRootGridIO == TRUE && ProblemType == 30)
	{
		if(PartitionNestedGrids)
		{
			printf("InitializeNew: Re-initialize NestedCosmologySimulation\n");
			if(NestedCosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL)
			{
				ENZO_FAIL("Error in NestedCosmologySimulationReInitialize.");
			}
		}
		else
		{
			printf("InitializeNew: Re-initialize CosmologySimulation\n");
			if(CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL)
			{
				ENZO_FAIL("Error in CosmologySimulationReInitialize.");
			}
		}
	}

	// For PhotonTest, using ParallelGridIO, initialize data only after
	// partitioning grid.  Last argument tells it it's the 2nd pass.

#ifdef TRANSFER
	if(ParallelRootGridIO == TRUE && ProblemType == 50)
		if(PhotonTestInitialize(fptr, Outfptr, TopGrid, MetaData, true) == FAIL)
			ENZO_FAIL("Error in PhotonTestInitialize(2nd pass).");
#endif

	if(ProblemType == 59)
		if(DrivenFlowInitialize(fptr, Outfptr, TopGrid, MetaData, 1) == FAIL)
			ENZO_FAIL("Error in DrivenFlowInitialize with SetBaryons");

	// For problem 60, using ParallelGridIO, read in data only after
	// partitioning grid.

	if(ParallelRootGridIO == TRUE && ProblemType == 60)
		if(TurbulenceSimulationReInitialize(&TopGrid, MetaData) == FAIL)
			ENZO_FAIL("Error in TurbulenceSimulationReInitialize.");

	if(ProblemType == 106)
	{
		if(TurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 1) == FAIL)
			ENZO_FAIL("Error in TurbulenceReInitialize.\n");
		//  if (HydroMethod == Zeus_Hydro) ConvertTotalEnergyToGasEnergy(&TopGrid);
	}

	if(ProblemType == 201)
		if(MHD2DTestInitialize(fptr, Outfptr, TopGrid, MetaData, 1) == FAIL)
			ENZO_FAIL("Error in MHD2DTestReInitialize.\n");

	if(ProblemType == 202)
		CollapseMHD3DInitialize(fptr, Outfptr, TopGrid, MetaData, 1);

	// For ProblemType 203 (MHD Turbulence Simulation we only initialize the data
	// once the topgrid has been split.
	if(ProblemType == 203)
		if(MHDTurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 1) == FAIL)
			ENZO_FAIL("Error in MHDTurbulenceReInitialize.\n");

	// initialize the data once the topgrid has been split.
	if(ProblemType == 210)
		if(MHDDecayingRandomFieldInitialize(fptr, Outfptr, TopGrid, MetaData, 1) == FAIL)
			ENZO_FAIL("Error in MHDDecayingRandomField ReInitialize.\n");

	if(ProblemType == 501)
		if(MHDProfileInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior, false) == FAIL)
			ENZO_FAIL("Error in MHDProfileInitialize with ParallelRootGridIO==true");

	PrintMemoryUsage("After 2nd pass");

	return SUCCESS;
}
