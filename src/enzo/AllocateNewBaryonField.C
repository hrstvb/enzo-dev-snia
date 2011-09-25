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
void *AllocateNewBaryonField(int size)
{
#ifndef MEMORY_POOL
  return new float[size];
#else
  return BaryonFieldMemoryPool->GetMemory(sizeof(float)*size);
#endif
}

void FreeBaryonFieldMemory(float *BF)
{
#ifndef MEMORY_POOL
  delete [] BF;
#else
  if (BF != NULL)
    BaryonFieldMemoryPool->FreeMemory(BF);
#endif
  return;
}

void FreeParticleMemory(void *PF)
{
#ifndef MEMORY_POOL
  delete [] PF;
#else
  return ParticleMemoryPool->FreeMemory(PF);
#endif
  return;
}

