/***********************************************************************
/
/  function AllocateNewBaryonField
/
/  written by: Tom Abel
/  date:       September 2011
/  modified:   
/
/  PURPOSE: Allocates space for the baryons. 
/           This can use the memory pool. 
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
float *AllocateNewBaryonField(int size)
{
#ifndef MEMORY_POOL
  return new float[size];
#else
  return static_cast<float*>(BaryonFieldMemoryPool->GetMemory(sizeof(float)*size));
#endif
}

void FreeBaryonFieldMemory(float *BF)
{
#ifndef MEMORY_POOL
    delete [] BF;
#else
    BaryonFieldMemoryPool->FreeMemory(BF);
#endif
  return;
}

int *AllocateNewFlaggingField(int size)
{
#ifndef MEMORY_POOL
  return new int[size];
#else
  //  fprintf(stderr,"Allocate new flagging field.");
  //  fflush(stderr);
  return static_cast<int*>(FlaggingFieldMemoryPool->GetMemory(sizeof(int)*size));
#endif
}

void FreeFlaggingFieldMemory(int *FF)
{
#ifndef MEMORY_POOL
  delete [] FF;
#else
  FlaggingFieldMemoryPool->FreeMemory(FF);
  //  fprintf(stderr, "Free flagging field.");
  //  fflush(stderr);
#endif
  return;
}


void FreeParticleMemory(void *PF)
{
#ifndef MEMORY_POOL
    delete [] PF;
#else
    ParticleMemoryPool->FreeMemory(PF);
#endif
  return;
}

