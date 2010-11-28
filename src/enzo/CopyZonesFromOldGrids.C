/***********************************************************************
/
/  COPY ZONES (non-blocking) FROM OLD TO NEW GRIDS (REBUILD HIERARCHY)
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
 
#include <stdio.h>
#include <string.h>
#include <time.h>
 
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "communication.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int CopyZonesFromOldGrids(LevelHierarchyEntry *OldGrids,
			  TopGridData *MetaData,
			  ChainingMeshStructure ChainingMesh)
{

  int i, grid1, NumberOfGrids, counter;
  FLOAT ZeroVector[] = {0,0,0};
  LevelHierarchyEntry *Temp;
  LevelHierarchyEntry **GridList;
  
  /* For threading, create a grid pointer array.  Can't use
     GenerateGridArray because GridHierarchyEntry isn't valid
     anymore. */

  NumberOfGrids = 0;
  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel)
    NumberOfGrids++;

  counter = 0;
  GridList = new LevelHierarchyEntry*[NumberOfGrids];
  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel)
    GridList[counter++] = Temp;

  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];

  /* Find sibling grids for all subgrids */

//#pragma omp parallel for schedule(static)
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    GridList[grid1]->GridData->FastSiblingLocatorFindSiblings
      (&ChainingMesh, &SiblingList[grid1],
       MetaData->LeftFaceBoundaryCondition,
       MetaData->RightFaceBoundaryCondition);

  /* Post receives, looping over the old grids */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

#pragma omp parallel for schedule(static) private(i)
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    // For each of the sibling grids, copy data.
    for (i = 0; i < SiblingList[grid1].NumberOfSiblings; i++)
      SiblingList[grid1].GridList[i]->CopyZonesFromGrid
	(GridList[grid1]->GridData, ZeroVector);

  }

  /* Send data */

  CommunicationDirection = COMMUNICATION_SEND;

#pragma omp parallel for schedule(static) private(i)
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    // For each of the sibling grids, copy data.
    for (i = 0; i < SiblingList[grid1].NumberOfSiblings; i++)
      SiblingList[grid1].GridList[i]->CopyZonesFromGrid
	(GridList[grid1]->GridData, ZeroVector);

    /* Delete all fields (only on the host processor -- we need
       BaryonField on the receiving processor) after sending them.  We
       only delete the grid object on all processors after
       everything's done. */

    if (GridList[grid1]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      GridList[grid1]->GridData->DeleteAllFields();

  }

  /* Receive data */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("");

  /* Delete old grids */

#pragma omp parallel for schedule(static)
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    delete [] SiblingList[grid1].GridList;
    delete GridList[grid1]->GridData;
    GridList[grid1]->GridData = NULL;
  }

  delete [] GridList;
  delete [] SiblingList;

  return SUCCESS;

}
