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
 
 
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    FreeParticleMemory(ParticleAcceleration[dim]);
    FreeBaryonFieldMemory(AccelerationField[dim]);
 
    ParticleAcceleration[dim]      = NULL;
    AccelerationField[dim]         = NULL;
  }
  //  FreeParticleMemory(ParticleAcceleration[MAX_DIMENSION]);
  //  ParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    FreeBaryonFieldMemory(OldBaryonField[i]);
    OldBaryonField[i] = NULL;
  }
 
  FreeBaryonFieldMemory(GravitatingMassField);
  FreeBaryonFieldMemory(GravitatingMassFieldParticles);
 
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;

#ifdef SAB
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    if (OldAccelerationField[dim] != NULL) {
      FreeBaryonFieldMemory(OldAccelerationField[dim]);
      OldAccelerationField[dim] = NULL;
    }
#endif

}
