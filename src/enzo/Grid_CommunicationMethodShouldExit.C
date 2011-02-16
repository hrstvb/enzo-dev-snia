/***********************************************************************
/
/  GRID CLASS: Communication helper function
/
/  written by: Greg Bryan
/  date:       
/  modified:   
/
/  PURPOSE: This is a simple helper function which determines if the 
/     method should return immediately because of communication-mode 
/     reasons.  Assumes that info is being sent from the "other-grid" 
/     processor to the "this-grid" processor (i.e. the one that holds 
/     this object). 
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "Parallel.h"

using Parallel::CommunicationDirection;

int grid::CommunicationMethodShouldExit(grid *OtherGrid) {

  /* Return if neither grid lives on this processor. */
  //    if (NumberOfProcessors == 1) return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != OtherGrid->ProcessorNumber)
    return SUCCESS;

  /* If the two grids are on the same processor then return if
     either in post-receive or receive modes to avoid duplicating method
     (i.e. action is only carried out if in send mode (or send-receive)). */

  if (ProcessorNumber == OtherGrid->ProcessorNumber)
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_RECEIVE)
      return SUCCESS;

  /* If in send mode then exit if the send grid is not on this processor. */

  if (CommunicationDirection == COMMUNICATION_SEND &&
      MyProcessorNumber != OtherGrid->ProcessorNumber)
    return SUCCESS;

  /* If in either receive phase then exit if receive grid is not on this 
     processor. */

  if ((CommunicationDirection == COMMUNICATION_RECEIVE ||
       CommunicationDirection == COMMUNICATION_POST_RECEIVE) &&
      MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  return FAIL; /* i.e. method should not exit immediately. */
}
