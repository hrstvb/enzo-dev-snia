/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This routine processes the receives stored in the 
/           CommunicationReceive stack.  Each receive is tagged with a 
/           type which indicates which method to call 
/           (and a record of the arguments).
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Parallel.h"
#include "CommunicationUtilities.h"

using namespace Parallel;
 
#ifdef USE_MPI
static MPI_Arg ListOfIndices[Parallel::MAX_RECEIVE_BUFFERS];
static MPI_Status ListOfStatuses[Parallel::MAX_RECEIVE_BUFFERS];
#endif /* USE_MPI */

double ReturnWallTime(void);

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[],
				int NumberOfSubgrids[],
				int FluxFlag, TopGridData* MetaData)
{

#ifdef USE_MPI

  int NoErrorSoFar = TRUE;
  int Zero[] = {0, 0, 0};
  FLOAT ZeroFL[] = {0, 0, 0};

  /* Set the communication mode. */

  Parallel::CommunicationDirection = COMMUNICATION_RECEIVE;

//    printf("P(%"ISYM") in CRH with %"ISYM" requests\n", MyProcessorNumber,
//	   CommunicationReceiveIndex);

  MPI_Arg NumberOfCompleteRequests, TotalReceives;
  int ReceivesCompletedToDate = 0, index, errcode, SUBling, level,
    igrid, isubgrid, dim, FromStart, FromNumber, ToStart, ToNumber;
  int GridDimension[MAX_DIMENSION];
  FLOAT EdgeOffset[MAX_DIMENSION];
  grid *grid_one, *grid_two;
  MPI_Request *AllRequests = NULL;
  MPIBuffer CurrentBuffer;
  GenerateMPIRequestArray(AllRequests);
  TotalReceives = CommunicationMPIBuffer.size();
  
#ifdef TRANSFER
  PhotonPackageEntry *PP;
#endif

  MPI_Arg error_code;
  char error_string[1024];
  MPI_Arg length_of_error_string, error_class;
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;

  while (ReceivesCompletedToDate < TotalReceives) {

    /* Call the MPI wait handler. */

    float time1 = ReturnWallTime();

//    printf("P%d ::: BEFORE MPI_Waitsome ::: %"ISYM" %"ISYM" %"ISYM"\n", 
//	   MyProcessorNumber,
//	   TotalReceives, ReceivesCompletedToDate, NumberOfCompleteRequests);
    MPI_Waitsome(TotalReceives, AllRequests,
		 &NumberOfCompleteRequests, ListOfIndices, ListOfStatuses);
//    printf("P%d: MPI: %"ISYM" %"ISYM" %"ISYM"\n", MyProcessorNumber,
//	   TotalReceives, 
//	   ReceivesCompletedToDate, NumberOfCompleteRequests);

    CommunicationTime += ReturnWallTime() - time1;

    /* Error check */

#if 0
    if (NumberOfCompleteRequests == MPI_UNDEFINED) {
      ENZO_FAIL("Error in MPI_Waitsome\n");
    }
#endif

    mpi_header CurrentHeader;

    /* Should loop over newly received completions and check error msgs now. */
    for (index = 0; index < NumberOfCompleteRequests; index++) {
      if (ListOfStatuses[index].MPI_ERROR != 0) {
	CurrentBuffer = GetMPIBuffer(ListOfIndices[index]);
	CurrentHeader = CurrentBuffer.ReturnHeader();
	if (NoErrorSoFar) {
	  fprintf(stderr, "MPI Error on processor %"ISYM". "
		  "Error number %"ISYM" on request %"ISYM"\n",
		  MyProcessorNumber, ListOfStatuses[index].MPI_ERROR, index);
	  NoErrorSoFar = FALSE;
	}
	fprintf(stdout, "P(%"ISYM") index %"ISYM" -- mpi error %"ISYM"\n", 
		MyProcessorNumber, ListOfIndices[index], 
		ListOfStatuses[index].MPI_ERROR);
	fprintf(stdout, "%"ISYM": Type = %d, Grid1 = %d, Request = %p\n",
		ListOfIndices[index], 
		CurrentHeader.CallType,
		CurrentBuffer.ReturnGridOne()->GetGridID(),
		AllRequests[ListOfIndices[index]]);
	fprintf(stdout, "%"ISYM": buffer = %p\n", ListOfIndices[index],
		CurrentBuffer);

	MPI_Error_class(ListOfStatuses[index].MPI_ERROR, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	fprintf(stderr, "P%d: %s\n", MyProcessorNumber, error_string);
	MPI_Error_string(ListOfStatuses[index].MPI_ERROR, error_string, &length_of_error_string);
	fprintf(stderr, "P%d: %s\n", MyProcessorNumber, error_string);
	ENZO_FAIL("");
      }
    } // ENDFOR index

    /* Loop over the receive handles, looking for completed (i.e. null)
       requests associated with unprocessed (i.e. non-null) grids. 
       It's insufficient to just loop over newly completed receives because
       there may be some completed receives which were not processed due
       to dependence issues. */

    for (index = 0; index < TotalReceives; index++) {

      CurrentBuffer = GetMPIBuffer(index);

      if (CurrentBuffer.ReturnGridOne() != NULL &&
	  AllRequests[index] == MPI_REQUEST_NULL) {

// 	if(CommunicationReceiveCallType[index]==2)
// 	  fprintf(stdout, "P%d: %d %d %d %d %d %d\n", MyProcessorNumber,
//		  index, 
// 		CommunicationReceiveCallType[index],
// 		CommunicationReceiveGridOne[index]->GetGridID(),
// 		CommunicationReceiveGridTwo[index]->GetGridID(),
// 		CommunicationReceiveMPI_Request[index],
// 		CommunicationReceiveDependsOn[index]);

	/* If this depends on an un-processed receive, then skip it. */
	/* JHW (Feb. 2011): Removed.  Should re-implement dependencies later. */

//	if (CommunicationReceiveDependsOn[index] != 
//	    COMMUNICATION_NO_DEPENDENCE)
//	  if (CommunicationReceiveGridOne[CommunicationReceiveDependsOn[index]]
//	      != NULL)
//	    continue;

	CurrentHeader = CurrentBuffer.ReturnHeader();
	grid_one = CurrentBuffer.ReturnGridOne();
	grid_two = CurrentBuffer.ReturnGridTwo();
	//CommunicationReceiveIndex = index;

	/* Copy out the argument if needed */

	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  EdgeOffset[dim] = CurrentHeader.FArg[dim];
	  
	/* Handle the buffers received, calling the appropriate method. */

	switch (CurrentHeader.CallType) {

	case 1:
	  errcode = grid_one->InterpolateBoundaryFromParent(grid_two, index);
	  break;

	case 2:
	  errcode = grid_one->CopyZonesFromGrid(grid_two, EdgeOffset, index);
	  break;

	case 3:
	  errcode = grid_one->
	    DepositParticlePositions(grid_two, CurrentHeader.FArg[0],
				     CurrentHeader.IArg[0], index);
	  break;

	case 4:
	  errcode = grid_one->CopyParentToGravitatingFieldBoundary
	    (grid_two, index);
	  break;

	case 5:
	  errcode = grid_one->DepositBaryons
	    (grid_two, CurrentHeader.FArg[0], index);
	  break;

	case 6:
	  errcode = grid_one->AddOverlappingParticleMassField
	    (grid_two, EdgeOffset, index);
	  break;

	case 7:
	  errcode = grid_one->PreparePotentialField(grid_two, index);
	  break;

	case 8:
	  errcode = grid_one->CopyOverlappingMassField(grid_two, EdgeOffset,
						       index);
	  break;

	case 9:
	  errcode = grid_one->CopyPotentialField(grid_two, EdgeOffset,
						 index);
	  break;

	case 10:
	  errcode = grid_one->InterpolateAccelerations(grid_two, index);
	  break;
      
	case 11:  /* Note this one involves two calls. */

	  /* Project subgrid's refined fluxes to the level of this grid. */

	  if (grid_one->GetProjectedBoundaryFluxes
	      (grid_two, 0, 0, SubgridFluxesRefined, index) == FAIL) {
	    ENZO_FAIL("Error in grid->GetProjectedBoundaryFluxes.\n");
	  }
	
	  /* Correct this grid for the refined fluxes (step #19)
	     (this also deletes the fields in SubgridFluxesRefined). */

	  // For SUBlings, the subgrid number is flagged by setting it to negative

	  igrid = CurrentHeader.IArg[0];
	  isubgrid = CurrentHeader.IArg[1];
	  SUBling = CurrentHeader.IArg[2];
#ifdef FLUX_FIX
	  if ((errcode = grid_two->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[igrid][isubgrid], &SubgridFluxesRefined, 
	       SubgridFluxesEstimate[igrid][NumberOfSubgrids[igrid] - 1],
	       SUBling, MetaData)) == FAIL) {
	    ENZO_FAIL("Error in grid->CorrectForRefinedFluxes.\n");
	  }
#else
	  if ((errcode = grid_two->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[igrid][isubgrid], &SubgridFluxesRefined, 
	       SubgridFluxesEstimate[igrid][NumberOfSubgrids[igrid] - 1]))
	      == FAIL) {
	    ENZO_FAIL("Error in grid->CorrectForRefinedFluxes.\n");
	  }
#endif
	  break;

	case 12:
	  errcode = grid_one->ProjectSolutionToParentGrid(*grid_two, index);
	  break;

	case 13:
	  errcode = grid_one->SetParticleMassFlaggingField(index);
	  break;

	case 14:
	  FromStart  = CurrentHeader.IArg[0];
	  FromNumber = CurrentHeader.IArg[1];
	  ToStart    = CurrentHeader.IArg[2];
	  errcode = grid_one->CommunicationSendParticles
	    (grid_two, MyProcessorNumber, FromStart, FromNumber, ToStart,
	     index);
	  break;

#ifdef TRANSFER
	case 15:
	  ToNumber = CurrentHeader.IArg[0];
	  FromNumber = CurrentHeader.IArg[1];
	  PP = grid_one->ReturnPhotonPackagePointer();
	  errcode = grid_one->CommunicationSendPhotonPackages
	    (grid_two, MyProcessorNumber, ToNumber, FromNumber,
	     &PP->NextPackage, index);
	  break;
#endif /* TRANSFER */

	case 16:
	  for (dim = 0; dim < MAX_DIMENSION; dim++)
	    GridDimension[dim] = CurrentHeader.IArg[dim];
	  // The last arguments (all zeros here) are only used in
	  // post-receive mode.
	  errcode = grid_one->CommunicationSendRegion
	    (grid_two, MyProcessorNumber, ALL_FIELDS, NEW_ONLY, Zero,
	     GridDimension, 0, NULL, NULL, ZeroFL, Zero, index);
	  break;

	case 17:
	  errcode = grid_one->InterpolateParticlesToGrid(NULL, index);
	  break;

	case 18:
	  errcode = grid_one->CommunicationSendStars
	    (grid_two, MyProcessorNumber, index);
	  break;

#ifdef TRANSFER
	case 19:
	  level = CurrentHeader.IArg[0];
	  errcode = grid_one->SetSubgridMarkerFromParent
	    (grid_two, level, index);
	  break;
#endif

	default:
	  ENZO_VFAIL("Unrecognized call type %"ISYM"\n",
		     CurrentHeader.CallType)

	} // end: switch on call type

	/* Report error if there has been one in any of the above calls. */

	if (errcode == FAIL) {
	  ENZO_VFAIL("Error in CommunicationReceiveHandler, method %"ISYM"\n",
		     CurrentHeader.CallType)

	}

	/* Mark this receive complete. */

	CurrentBuffer.MarkComplete();
	//	MPI_Request_free(CommunicationReceiveMPI_Request+index);
	ReceivesCompletedToDate++;

      } // end: if statement to check if receive should be processed

    } // end: loop over all receives

  } // end: while loop waiting for all receives to be processed

  Parallel::CommunicationReceiveIndex = 0;

  /* Reset the communication mode. */

  Parallel::CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  //  printf("P(%d) out of CRH\n", MyProcessorNumber);

#endif /* USE_MPI */

  return SUCCESS;

}
