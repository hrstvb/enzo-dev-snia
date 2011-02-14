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
#include "Parallel.h"

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
namespace Parallel {

  int CreateMPIHeaderType(MPI_Datatype &HeaderType)
  {

    int i;
    mpi_header header;
    MPI_Datatype type[6] = {MPI_CHAR, IntDataType, IntDataType, IntDataType,
			    FLOATDataType, IntDataType};
    int blocklen[6] = {1, 2, 3, 3, 3, 3};
    MPI_Aint disp[6];
    MPI_Aint base;
    
    /* Compute the displacements inside the struct, then create header
       type */

    MPI_Address(&header, disp);
    MPI_Address(header.GridNum, disp+1);
    MPI_Address(header.GridOffset, disp+2);
    MPI_Address(header.GridDims, disp+3);
    MPI_Address(header.FArg, disp+4);
    MPI_Address(header.IArg, disp+5);
    base = disp[0];
    for (i = 0; i < 6; i++) disp[i] -= base;
    MPI_Type_struct(6, blocklen, disp, type, &HeaderType);
    MPI_Type_commit(&HeaderType);

    CommunicationHeaderType = HeaderType;

    return SUCCESS;

  }

  /*----------------------------------------------------------------------*/

  MPIBuffer::MPIBuffer(const int grid1, 
		       const int grid2, 
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

    request = MPI_REQUEST_NULL;
    tag = MPI_Tag;

    /* Create header */
  
    mpi_header header;
    header.CallType = CallType;
    header.GridNum[0] = grid1;
    header.GridNum[1] = grid2;

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

    MPI_Datatype type[2] = {CommunicationHeaderType, BufferType};
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
    MPI_Recv(this->message, 1, MessageType, FromProcessor, this->tag,
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return SUCCESS;
  }

  int MPIBuffer::IRecvBuffer(int FromProcessor)
  {
    MPI_Irecv(this->message, 1, MessageType, FromProcessor, this->tag,
	      MPI_COMM_WORLD, &this->request);
    return SUCCESS;
  }

  int MPIBuffer::SendBuffer(int ToProcessor)
  {
    CommunicationBufferedSend(this->message, 1, this->MessageType,
			      ToProcessor, this->tag, MPI_COMM_WORLD,
			      BUFFER_IN_PLACE);
    return SUCCESS;
  }

} // END NAMESPACE
#endif /* USE_MPI */
