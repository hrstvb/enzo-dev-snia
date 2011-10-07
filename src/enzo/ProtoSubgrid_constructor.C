/***********************************************************************
/
/  PROTOSUBGRID CLASS (CONSTRUCTOR)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void FreeFlaggingFieldMemory(int *FF);

 
ProtoSubgrid::ProtoSubgrid()
{
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    GridDimension[dim] = 1;
    StartIndex[dim]    = 0;
    EndIndex[dim]      = 0;
    Signature[dim]     = NULL;
  }
 
  GridFlaggingField = NULL;
 
  NumberFlagged = INT_UNDEFINED;
}
 
#ifdef MEMORY_POOL
void* ProtoSubgrid::operator new(size_t object_size)
{
  return ProtoSubgridMemoryPool->GetMemory(object_size);
}

void ProtoSubgrid::operator delete(void* object)
{
  ProtoSubgridMemoryPool->FreeMemory(object);
  return;
}
#endif
 
 
ProtoSubgrid::~ProtoSubgrid()
{
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    {
      FreeFlaggingFieldMemory(Signature[dim]);
      Signature[dim] = NULL;
    }
  //    delete [] Signature[dim];
  {
  FreeFlaggingFieldMemory(GridFlaggingField);
  GridFlaggingField = NULL;
  }
    // delete [] GridFlaggingField;
}

