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
#include "Grid.h"
#include "CommunicationUtilities.h"

int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);

int InitializeNew_PartitionGrids(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior, float *Initialdt, FILE *fptr, FILE *Outfptr, bool partitionNestedGrids)
{
	HierarchyEntry *g = &TopGrid;
	int level = 0;

	if(partitionNestedGrids)
		TRACEF("Partitioning all level(s) ...");
	else
		TRACEF("Partitioning TopGrid (level 0) only ...");

	while(g)
	{
		if(g->NextGridThisLevel)
		{
			TRACEF("  Level %" ISYM " already partitioned.", level);
		}
		else
		{
			TRACEF("  Partitioning single initial grid on level %" ISYM " ...", level);
			CommunicationPartitionGrid(g, level);
		}

		if(!partitionNestedGrids)
			break;

		level++;
		g = g->NextGridNextLevel;
	}

	TRACEF("Partitioning done.");
	return SUCCESS;
}
