/***********************************************************************
/
/  GRID CLASS (REMOVE ALL FIELDS)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
void FreeBaryonFieldMemory(float *BF);

void grid::DeleteAllFields()
{
 
  int i;
 
  this->DeleteParticles();
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    FreeParticleMemory(ParticleAcceleration[i]);
    FreeBaryonFieldMemory(AccelerationField[i]);
 
    ParticleAcceleration[i]      = NULL;
    AccelerationField[i]         = NULL;
  }
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    FreeBaryonFieldMemory(BaryonField[i]);
    FreeBaryonFieldMemory(OldBaryonField[i]);
    BaryonField[i]    = NULL;
    OldBaryonField[i] = NULL;
  }

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      FreeBaryonFieldMemory(OldAccelerationField[i]);
      OldAccelerationField[i] = NULL;
    }
#endif
 
  FreeBaryonFieldMemory(PotentialField);
  FreeBaryonFieldMemory(GravitatingMassField);
  delete [] GravitatingMassFieldParticles;
 
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
 
}
