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

#include "DebugMacros.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "Grid.h"

int InitializeNew_TopGrid(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr)
{
	// Error check the rank
	//printf("This should only run if not a restart!");
	if(MetaData.TopGridRank < 0 || MetaData.TopGridRank > 3)
	{
		ENZO_VFAIL("TopGridRank = %"ISYM" ill defined.\n", MetaData.TopGridRank)
	}

	// Error check the dimensions and at the same time add ghost zones
	for(int dim = 0; dim < MetaData.TopGridRank; dim++)
	{
		if(MetaData.TopGridDims[dim] < 1 || MetaData.TopGridDims[dim] > 8192)
		{
			ENZO_VFAIL("TopGridDims[%"ISYM"] = %"ISYM" ill defined.\n", dim, MetaData.TopGridDims[dim])
		}
		MetaData.TopGridDims[dim] =
				(MetaData.TopGridDims[dim] > 1) ? MetaData.TopGridDims[dim] + 2 * NumberOfGhostZones : 1;
	}

	// Create the top grid, prepare it, set the time and parameters
	TopGrid.GridData = new grid;

	TopGrid.GridData->PrepareGrid(MetaData.TopGridRank, MetaData.TopGridDims, DomainLeftEdge, DomainRightEdge,
									MetaData.NumberOfParticles);
	TopGrid.GridData->SetTime(MetaData.Time);
	TopGrid.GridData->SetHydroParameters(MetaData.CourantSafetyNumber, MetaData.PPMFlatteningParameter,
											MetaData.PPMDiffusionParameter, MetaData.PPMSteepeningParameter);
	TopGrid.GridData->SetGravityParameters(MetaData.GravityBoundary);

	// Repair TopGridDims (subtract ghost zones added earlier)
	for(int dim = 0; dim < MetaData.TopGridRank; dim++)
		MetaData.TopGridDims[dim] = max(MetaData.TopGridDims[dim] - 2 * NumberOfGhostZones, 1);

	// Set TopGrid Hierarchy Entry
	TopGrid.NextGridThisLevel = NULL;  // always true
	TopGrid.ParentGrid = NULL;  // always true
	TopGrid.NextGridNextLevel = NULL;  // can be reset by initializer

	return SUCCESS;
}
