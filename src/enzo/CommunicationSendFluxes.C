/***********************************************************************
/
/  COMMUNICATION ROUTINE: SEND FLUXES TO ANOTHER PROCESSOR
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
 
#include "ErrorExceptions.h"
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
#include "Parallel.h"

using namespace Parallel;
 
/* function prototypes */
 
int CommunicationSendFluxes(MPIBuffer *mbuffer, fluxes *Fluxes, int ToProc, 
			    int NumberOfFields, int Rank)
{
 
  /* Count space and Allocate buffer. */
 
  int dim1, dim2, field, i, TotalSize = 0, Sizes[MAX_DIMENSION], TempDim;
  for (dim1 = 0; dim1 < Rank; dim1++) {
    int size = 1;
    for (dim2 = 0; dim2 < Rank; dim2++) {
      TempDim = (Fluxes->LeftFluxEndGlobalIndex[dim1][dim2] -
	         Fluxes->LeftFluxStartGlobalIndex[dim1][dim2]) + 1;
      if (dim2 == dim1)
	TempDim = 1;
      size *= TempDim;
    }
    Sizes[dim1] = size;
    TotalSize += 2*size;
  }
 
  TotalSize *= NumberOfFields;
  float *buffer = new float[TotalSize];
 
  /* Pack buffer. */
 
  int index = 0;
  for (dim1 = 0; dim1 < Rank; dim1++)
    for (field = 0; field < NumberOfFields; field++) {
      for (i = 0; i < Sizes[dim1]; i++)
	buffer[index++] = Fluxes->LeftFluxes[field][dim1][i];
      for (i = 0; i < Sizes[dim1]; i++)
	buffer[index++] = Fluxes->RightFluxes[field][dim1][i];
    }
 
  /* send. */
 
#ifdef USE_MPI

#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  mbuffer->FillBuffer(FloatDataType, TotalSize, buffer);
  mbuffer->SendBuffer(ToProc);
 
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  endtime = MPI_Wtime();
  timer[12] += endtime-starttime;
  counter[12] ++;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
//  delete [] buffer;
 
  return SUCCESS;
}
