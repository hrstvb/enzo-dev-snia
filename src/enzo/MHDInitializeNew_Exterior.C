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
#include "Grid.h"

int InitializeNew_Exterior(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior,
float *Initialdt, FILE *fptr, FILE *Outfptr)
{
	// Initialize the exterior (unless it was set in the problem initializer)
	if(Exterior.AmIPrepared() != FALSE)
	{
		if(debug)
			printf("Exterior already initialized.\n");
		return SUCCESS;
	}

	if(debug)
		printf("Initializing Exterior ...\n");

	Exterior.Prepare(TopGrid.GridData);   // set rank and dims

	if(MetaData.BoundaryConditionName != NULL)
	{
		FILE *BCfptr;
		if((BCfptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL)
			ENZO_VFAIL("Error opening BC file: %s\n", MetaData.BoundaryConditionName)

		fprintf(stderr, "Opened BC file mode r\n");

		if(Exterior.ReadExternalBoundary(BCfptr) == FAIL)
			ENZO_FAIL("Error in ReadExternalBoundary.");

		fclose(BCfptr);
	}
	else
	{
		if(debug)
			fprintf(stderr, "InitializeExternalBoundaryFace\n");

		SimpleConstantBoundary = TRUE;
		for(int dim = 0; dim < MetaData.TopGridRank; dim++)
		{
			if(MetaData.LeftFaceBoundaryCondition[dim] != periodic
					|| MetaData.RightFaceBoundaryCondition[dim] != periodic)
			{
				SimpleConstantBoundary = FALSE;
			}
		}

		if(debug)
			fprintf(stderr, "SimpleConstantBoundary = %ISYM\n", SimpleConstantBoundary);

		float Dummy[TopGrid.GridData->ReturnNumberOfBaryonFields()];
		for(int fieldIndex = 0; fieldIndex < TopGrid.GridData->ReturnNumberOfBaryonFields(); fieldIndex++)
			Dummy[fieldIndex] = 0.0;

		for(int dim = 0; dim < MetaData.TopGridRank; dim++)
		{
			if(Exterior.InitializeExternalBoundaryFace(dim, MetaData.LeftFaceBoundaryCondition[dim],
														MetaData.RightFaceBoundaryCondition[dim], Dummy, Dummy) == FAIL)
			{
				ENZO_FAIL("Error in InitializeExternalBoundaryFace.");
			}
		}

		// Initialize particle boundary conditions
		Exterior.InitializeExternalBoundaryParticles(MetaData.ParticleBoundaryType);
	}  // end: if (MetaData.BoundaryConditionName != NULL)

	if(debug)
		printf("Initializing Exterior done.\n");

	return SUCCESS;
}
