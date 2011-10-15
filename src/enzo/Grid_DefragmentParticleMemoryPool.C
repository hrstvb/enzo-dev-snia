/***********************************************************************
/
/  GRID CLASS (DEFRAGMENTS PARTICLE MEMORY POOL)
/
/  written by: Tom Abel
/  date:       October 2011
/  modified1:
/
/  PURPOSE:
/       Copy all particle fields to a new baryon pool getting rid of wasted 
/       space.
/       Should only be called from DefragmentMemoryPools.C
/ 
************************************************************************/
 
//  Compute derived quantites
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::DefragmentParticleMemoryPool()
{
#ifdef MEMORY_POOL
  /*  Particles only: set up field quantities and allocate fields
      (we assume here the grid is uniform in each dimension) */
 

      // use Particle Memory Pool to allocate memory
  if (ParticleMass != NULL)
    {
      float *temp = static_cast<float*>(ParticleMemoryPool->GetMemory(sizeof(float)*NumberOfParticles));
      for (int index = 0; index<NumberOfParticles; index++ )
	temp[index] = ParticleMass[index];
      FreeParticleMemory(ParticleMass);
      ParticleMass = temp;
    }
	   
  if (ParticleNumber != NULL)
    {
      PINT *temp = static_cast<PINT*>(ParticleMemoryPool->GetMemory(sizeof(PINT)*NumberOfParticles));
      for (int index = 0; index<NumberOfParticles; index++ )
	temp[index] = ParticleNumber[index];
      FreeParticleMemory(ParticleNumber);
      ParticleNumber = temp;
    }

  if (ParticleType != NULL)
    {
      int *temp = static_cast<int*>(ParticleMemoryPool->GetMemory(sizeof(int)*NumberOfParticles));
      for (int index = 0; index<NumberOfParticles; index++ )
	temp[index] = ParticleType[index];
      FreeParticleMemory(ParticleType);
      ParticleType = temp;
    }

  for (int dim = 0; dim < GridRank; dim++) {
    if (ParticlePosition[dim] != NULL)
      {
	FLOAT *temp = static_cast<FLOAT*>(ParticleMemoryPool->GetMemory(sizeof(FLOAT)*NumberOfParticles));
	for (int index = 0; index<NumberOfParticles; index++ )
	  temp[index] = ParticlePosition[dim][index];
	FreeParticleMemory(ParticlePosition[dim]);
	ParticlePosition[dim] = temp;
      }
    if (ParticleVelocity[dim] != NULL)
      {
	float *temp = static_cast<float*>(ParticleMemoryPool->GetMemory(sizeof(float)*NumberOfParticles));
	for (int index = 0; index<NumberOfParticles; index++ )
	  temp[index] = ParticleVelocity[dim][index];
	FreeParticleMemory(ParticleVelocity[dim]);
	ParticleVelocity[dim] = temp;
      }
  }
  for (int i = 0; i < NumberOfParticleAttributes; i++)
    if (ParticleAttribute[i] != NULL)
      {
	float *temp = static_cast<float*>(ParticleMemoryPool->GetMemory(sizeof(float)*NumberOfParticles));
	for (int index = 0; index<NumberOfParticles; index++ )
	  temp[index] = ParticleAttribute[i][index];
	FreeParticleMemory(ParticleAttribute[i]);
	ParticleAttribute[i] = temp;
      }
  
#endif // MEMORY_POOL  
  return;
}





 

 
