/***********************************************************************
/
/  GRID CLASS (BEFORE REBUILDING, REMOVED UNNEEDED ARRAYS)
/
/  written by: Greg Bryan
/  date:       June, 1995
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
 
void grid::CleanUp()
{
 
  int i;
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    FreeParticleMemory(ParticleAcceleration[i]);
//    delete [] AccelerationField[i];
 
    ParticleAcceleration[i]      = NULL;
//    AccelerationField[i]         = NULL;
  }
  FreeParticleMemory(ParticleAcceleration[MAX_DIMENSION]);
  ParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    FreeBaryonFieldMemory(OldBaryonField[i]);
    OldBaryonField[i] = NULL;
  }
 
  FreeBaryonFieldMemory(GravitatingMassField);
  FreeBaryonFieldMemory(GravitatingMassFieldParticles);
 
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      FreeBaryonFieldMemory(OldAccelerationField[i]);
      OldAccelerationField[i] = NULL;
    }
#endif

}
