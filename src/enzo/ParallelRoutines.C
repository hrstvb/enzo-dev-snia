/***********************************************************************
/
/  PARALLEL ROUTINES
/
/  written by: John Wise
/  date:       February, 2011
/  modified1:  
/
***********************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "GroupPhotonList.h"
#include "Parallel.h"

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
namespace Parallel {

  /*----------------------------------------------------------------------*/

  void GenerateMPIRequestArray(MPI_Request *result)
  {
    
    int i;
    int size = CommunicationMPIBuffer.size();
    result = new MPI_Request[size];
    list<MPIBuffer>::iterator it;
    for (it = CommunicationMPIBuffer.begin();
	 it != CommunicationMPIBuffer.end(); ++it) {
      i = (*it).ReturnIndex();
      result[i] = (*it).ReturnRequest();
    }

  }

  /*----------------------------------------------------------------------*/

  MPIBuffer GetMPIBuffer(int num)
  {
    list<MPIBuffer>::iterator it = CommunicationMPIBuffer.begin();
    advance(it,num);
    return *it;
  }

  /*----------------------------------------------------------------------*/

  MPIBuffer::MPIBuffer(void)
  {
    this->BufferType = MPI_DATATYPE_NULL;
    this->MessageType = MPI_DATATYPE_NULL;
    this->request = MPI_REQUEST_NULL;
    this->tag = 0;
    this->grid_one = NULL;
    this->grid_two = NULL;
    this->message = new enzo_message;
  }

  MPIBuffer::MPIBuffer(grid *grid1, 
		       grid *grid2, 
		       const int CallType,
		       const int MPI_Tag,
		       const int GridOffset[3],
		       const int GridDims[3],
		       const FLOAT FloatArgs[3],
		       const int IntArgs[3])
  {

    int i;

#pragma omp critical
    {
      this->message_index = CommunicationReceiveIndex++;
    }

    this->request = MPI_REQUEST_NULL;
    this->tag = MPI_Tag;
    this->grid_one = grid1;
    this->grid_two = grid2;

    /* Create header */
  
    mpi_header header;
    header.CallType = CallType;
    if (grid1 != NULL)
      header.GridNum[0] = grid1->GetGridID();
    else
      header.GridNum[0] = 0;
    if (grid2 != NULL)
      header.GridNum[1] = grid2->GetGridID();
    else
      header.GridNum[1] = 0;

    if (GridOffset != NULL) 
      for (i = 0; i < 3; i++) header.GridOffset[i] = GridOffset[i];
    else
      for (i = 0; i < 3; i++) header.GridOffset[i] = 0;

    if (GridDims != NULL) 
      for (i = 0; i < 3; i++) header.GridDims[i] = GridDims[i];
    else
      for (i = 0; i < 3; i++) header.GridDims[i] = 1;

    if (FloatArgs != NULL) 
      for (i = 0; i < 3; i++) header.FArg[i] = FloatArgs[i];
    else
      for (i = 0; i < 3; i++) header.FArg[i] = 0.0;

    if (IntArgs != NULL) 
      for (i = 0; i < 3; i++) header.IArg[i] = IntArgs[i];
    else
      for (i = 0; i < 3; i++) header.IArg[i] = 0;

    this->message = new enzo_message;
    this->message->header = header;

  }

  MPIBuffer::~MPIBuffer(void)
  {
    MPI_Type_free(&BufferType);
    MPI_Type_free(&MessageType);
    //delete [] this->message.buffer;
  }

  int MPIBuffer::FillBuffer(const MPI_Datatype BufferDataType, 
			    const int BufferSize,
			    void *buffer)
  {

    /* Create buffer as a contiguous datatype.  Note: it may be useful
       if we use vector datatypes in the future to provide more
       flexibility and less copying. */
  
    MPI_Type_contiguous(BufferSize, BufferDataType, &BufferType);
    MPI_Type_commit(&BufferType);

    /* Combine header and buffer, then register datatype */

    this->message->buffer = buffer;

    MPI_Datatype type[2] = {MPI_Header, BufferType};
    int blocklen[2] = {1,1};
    MPI_Aint disp[2];
    MPI_Aint base;

    MPI_Address(&(this->message), disp);
    MPI_Address(this->message->buffer, disp+1);
    base = disp[0];
    for (int i = 0; i < 2; i++) disp[i] -= base;
    MPI_Type_struct(2, blocklen, disp, type, &MessageType);
    MPI_Type_commit(&MessageType);

    return SUCCESS;

  }  

  int MPIBuffer::RecvBuffer(int FromProcessor)
  {
    MPI_Recv(this->message, 1, MessageType, (MPI_Arg) FromProcessor, 
	     this->tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return SUCCESS;
  }

  int MPIBuffer::IRecvBuffer(int FromProcessor)
  {
    MPI_Irecv(this->message, 1, MessageType, (MPI_Arg) FromProcessor, 
	      this->tag, MPI_COMM_WORLD, &this->request);
    // Append to receive buffer list
#pragma omp critical
    {
      CommunicationMPIBuffer.push_front(*this);
    }
    return SUCCESS;
  }

  int MPIBuffer::SendBuffer(int ToProcessor, int size)
  {
    CommunicationBufferedSend(this->message, 1, this->MessageType,
			      ToProcessor, this->tag, MPI_COMM_WORLD,
			      size);
    return SUCCESS;
  }

} // END NAMESPACE
#endif /* USE_MPI */
