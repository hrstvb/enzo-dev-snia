/***********************************************************************
/
/  COMMUNICATION SPECIFC GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       August, 2003
/  modified1:  February, 2011 (JHW) turning into a namespace
/
/  PURPOSE:
/    This is data required for the optimised communication routines.
/    It holds temporary data associated mostly with receive buffers.
/
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/
#ifndef __PARALLEL_H
#define __PARALLEL_H

#define COMMUNICATION_NO_DEPENDENCE -1

#include <stdlib.h>
#include <list>

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include "mpi_typedef.h"

using namespace std;

namespace Parallel
{

/* Set maximum number of receive buffers */

  static const size_t MAX_RECEIVE_BUFFERS = 50000;

/* Set the code for the no dependence for the DependsOn element (see below). */

/* This is the current mode for communiction.  There are two different ways
   that communication can be done.  The first is generally slower and is
   included only for error-checking; it is used by the routines in
   EvolveLevelRoutines.C which is generally not included in the Makefile
   (in favour of EvolveLevelRoutinesOptimiezed.C).

   1) single-phase: COMMUNICATION_SEND_RECEIVE
      All communication is carried out immediately.  That is, each 
      communication method is only called once and the function waits until
      the data actually arrives.

   2) triple-phase:  In this technique, the call for receiving the data comes
      first (COMMUNICATION_POST_RECEIVE) without waiting for the data to
      arrive.  This generates a buffer into which the data will go and it
      generates a CommunicationReceive handler, as detailed below.

      Then, the methods are all called again, with CommunicationDirection now
      set to COMMUNICATION_SEND, which generates the actual (buffered) send
      calls.  Once again, the call completes immediately and the routine
      CommunicationBufferedSend handles deallocating the buffers when the
      send completes.

      In the final phase, MPI_WAITSOME is used to find out receive buffers
      which have completed.  This then generates a third call to the
      appropriate method, this time with CommunicationDirection set to
      COMMUNICATION_RECEIVE.  There is a generalized routine
      (CommunicationReceiveHandler) which does this for all the receieve
      methods. */

  int CommunicationDirection;

/* This variable contains the most recent receive dependence; that is, the
   index of the receive handler which must complete first. */

  int CommunicationReceiveCurrentDependsOn;

/* This is the index of the current receive buffer.  Alterations of
   CommunicationReceiveIndex or CommunicationGridID MUST be done in
   omp critical sections because they are global. */

  int CommunicationReceiveIndex;

/* This is the grid ID of the receiving and sending grid that is being
   processed for use with thread-safe MPI tags.  Used in MPI
   headers. */

  //int CommunicationGridID[2];

/* More flags to ensure unique MPI tags where the receiving and
   sending grids are the same.  Necessary when copying overlapping
   regions. */

  //int CommunicationTags[3];

/* The following variables contain information about each receive buffer
   handler.  They are:

   CallType - an integer representing the method type that generated the
              receive handler in the first place and must be called to
	      complete the receive.
   MPI_Request - this is a pointer to the MPI receive request handle.
   GridOne - This is a pointer to the grid object which generated the handle.
   GridTwo - This is a pointer to the grid object which is an argument to
             the method.
   DependsOn - This is an index to the receive handle which must complete
               before this call can be processed.
   CallArgument - This is the value to any extra argument in the grid
                  method which generated the handle.                     */

#ifdef USE_MPI
  MPI_Datatype CommunicationHeaderType;
  //MPI_Request  CommunicationReceiveMPI_Request[MAX_RECEIVE_BUFFERS];
  //MPI_Datatype CommunicationReceiveMPI_Datatype[MAX_RECEIVE_BUFFERS];
  //void       *CommunicationReceiveBuffer[MAX_RECEIVE_BUFFERS];
  //int          CommunicationReceiveCallType[MAX_RECEIVE_BUFFERS];
  //grid        *CommunicationReceiveGridOne[MAX_RECEIVE_BUFFERS];
  //grid        *CommunicationReceiveGridTwo[MAX_RECEIVE_BUFFERS];
  //int          CommunicationReceiveDependsOn[MAX_RECEIVE_BUFFERS];
  //FLOAT CommunicationReceiveArgument[MAX_DIMENSION][MAX_RECEIVE_BUFFERS];
  //int CommunicationReceiveArgumentInt[MAX_DIMENSION][MAX_RECEIVE_BUFFERS];

  int CreateMPIHeaderType(MPI_Datatype &HeaderType);

  class MPIBuffer {

  private:
    int message_index;
    MPI_Datatype MessageType;
    MPI_Datatype BufferType;
    MPI_Request request;
    MPI_Arg tag;
    enzo_message *message;

  public:

    MPIBuffer(const int grid1, const int grid2, const int CallType,
	      const int MPI_Tag,
	      const int GridOffset[] = NULL,
	      const int GridDims[] = NULL,
	      const FLOAT FloatArgs[] = NULL,
	      const int IntArgs[] = NULL);

    ~MPIBuffer(void);
    
    int RecvBuffer(int FromProcessor);
    int IRecvBuffer(int FromProcessor);
    int SendBuffer(int ToProcessor);
    int FillBuffer(const MPI_Datatype BufferDataType, 
		   const int BufferSize,
		   void *buffer);
    
  }; // ENDCLASS    

  list<MPIBuffer> CommunicationMPIBuffer;

#endif /* USE_MPI */

} // END NAMESPACE
#endif  /* ifndef __PARALLEL_H */
