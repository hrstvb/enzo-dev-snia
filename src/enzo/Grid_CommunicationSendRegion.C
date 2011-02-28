/***********************************************************************
/
/  GRID CLASS (SEND FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#ifdef _OPENMP
#include "omp.h"
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Parallel.h"
#include "CommunicationUtilities.h"

using namespace Parallel;

// function prototypes
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */

 
int grid::CommunicationSendRegion(grid *ToGrid, int ToProcessor,int SendField,
				  int NewOrOld, int RegionStart[], int RegionDim[],
				  int CommType, grid* grid_one, grid* grid_two,
				  FLOAT CommArg[], int CommArgInt[],
				  int CommunicationIndex)
{
#ifdef USE_MPI 

  /* Return if not on processor. */

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

//  if (MyProcessorNumber != ProcessorNumber && 
//      MyProcessorNumber != ToProcessor)
//    return SUCCESS;

  int index, field, dim, Zero[] = {0, 0, 0};
 
  // Compute size of region to transfer
 
  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
 
  if (SendField == ACCELERATION_FIELDS)
    NumberOfFields = GridRank;
 
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
  int TransferSize = RegionSize * NumberOfFields;

  // Allocate MPI buffer
  MPIBuffer *mbuffer = NULL;
  if (Parallel::CommunicationDirection != COMMUNICATION_RECEIVE)
    mbuffer = new MPIBuffer(grid_one, grid_two, CommType, 
			    MPI_SENDREGION_TAG, RegionStart, RegionDim,
			    CommArg, CommArgInt);
  else {
    mbuffer = GetMPIBuffer(CommunicationIndex);  // Grab from list.
  }


  // Allocate buffer
  float *buffer = NULL;
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (float*) mbuffer->ReturnBuffer();
  else	   
    buffer = new float[TransferSize];
 
  // If this is the from processor, pack fields
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    //printf("SendRegion: RegionStart = %"ISYM" %"ISYM" %"ISYM"\n", RegionStart[0], RegionStart[1], RegionStart[2]);
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(BaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(OldBaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
      }
 
    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES)
      FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, buffer,
			   GravitatingMassFieldParticlesDimension,
			   GravitatingMassFieldParticlesDimension+1,
			   GravitatingMassFieldParticlesDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == GRAVITATING_MASS_FIELD)
      FORTRAN_NAME(copy3d)(GravitatingMassField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == POTENTIAL_FIELD)
      FORTRAN_NAME(copy3d)(PotentialField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	FORTRAN_NAME(copy3d)(AccelerationField[dim], &buffer[index],
			     GridDimension, GridDimension+1, GridDimension+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     RegionStart, RegionStart+1, RegionStart+2);
	index += RegionSize;
      }
  }

  /* Send buffer */
 
  /* Only send if processor numbers are not identical */
 
  if (ProcessorNumber != ToProcessor) {

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif

    /* Send the data if on send processor, but leave buffer until the data
       has been transfered out. */

    if (CommunicationDirection != COMMUNICATION_RECEIVE)
      mbuffer->FillBuffer(FloatDataType, TransferSize, buffer);

    if (MyProcessorNumber == ProcessorNumber) {
#ifdef MPI_INSTRUMENTATION
      if (traceMPI) 
	fprintf(tracePtr, "CSR Sending %"ISYM" floats from %"ISYM" to %"ISYM"\n", 
		TransferSize, MyProcessorNumber, ToProcessor);
#endif
      printf("CSR Sending %d bytes from %"ISYM" to %"ISYM"\n", 
	     TransferSize*sizeof(float), MyProcessorNumber, ToProcessor);
      mbuffer->SendBuffer(ToProcessor);
    }

    if (MyProcessorNumber == ToProcessor) {

      fprintf(stderr, "Waiting for %d floats at %d from %d\n", TransferSize, 
	      MyProcessorNumber, ProcessorNumber);

      /* Post the receive message without waiting for the message to
	 be received.  When the data arrives, this will be called again
	 in (the real) receive mode. */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {

	printf("Posting receive from P%"ISYM" for %"ISYM" floats in "
	       "comm index %"ISYM"\n", ProcessorNumber, 
	       TransferSize, mbuffer->ReturnIndex());
	mbuffer->IRecvBuffer(ProcessorNumber);
      }

      /* If in send-receive mode, then wait for the message now. */

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	mbuffer->RecvBuffer(ProcessorNumber);

    } // ENDIF ToProcessor

#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[5] += endtime-starttime;
    counter[5] ++;
    timer[6] += double(TransferSize);
    RecvComm += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
  } // ENDIF different processors
 
  /* If this is the to processor, and we're either in send-receive mode
     or receive mode, then unpack the data. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    if (ToProcessor != ProcessorNumber)
      fprintf(stderr, "Received %d floats at %d from %d\n", TransferSize, 
	      MyProcessorNumber, ProcessorNumber);

    index = 0;

    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->BaryonField[field];
	  ToGrid->BaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->OldBaryonField[field];
	  ToGrid->OldBaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}
 
    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
      delete ToGrid->GravitatingMassFieldParticles;
      ToGrid->GravitatingMassFieldParticles = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassFieldParticles,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == GRAVITATING_MASS_FIELD) {
      delete ToGrid->GravitatingMassField;
      ToGrid->GravitatingMassField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == POTENTIAL_FIELD) {
      delete ToGrid->PotentialField;
      ToGrid->PotentialField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->PotentialField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	delete ToGrid->AccelerationField[dim];
	ToGrid->AccelerationField[dim] = new float[RegionSize];
	FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->AccelerationField[dim],
			     RegionDim, RegionDim+1, RegionDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     Zero, Zero+1, Zero+2);
	index += RegionSize;
      }

    /* Only delete the buffer if we're in receive mode (in send mode
       it will be deleted by CommunicationBufferedSend and if we're in
       post-receive mode then it will be deleted when we get to
       receive-mode). */

    delete [] buffer;
			  
  } // ENDIF unpack
 
#endif /* USE_MPI */ 

  return SUCCESS;
}
