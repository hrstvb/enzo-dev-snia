/***********************************************************************
/
/  GRID CLASS (DESTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//
//  Grid destructor
//
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
 
void DeleteFluxes(fluxes *Fluxes);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
void DeleteStarList(Star * &Node);
#ifdef TRANSFER
PhotonPackageEntry*  DeletePhotonPackage(PhotonPackageEntry *PP);
#endif /* TRANSFER */
 
grid::~grid()
{
 
  int i;
 
  /* Error check. */
 
#ifdef UNUSED
  if (NumberOfParticles > 0) {
    fprintf(stderr, "warning: destroying live particles (%"ISYM").\n",
	    NumberOfParticles);
  /* exit(EXIT_FAILURE); */
  }
#endif /* UNUSED */
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] CellLeftEdge[i];
    delete [] CellWidth[i];
    FreeParticleMemory(ParticlePosition[i]);
    FreeParticleMemory(ParticleVelocity[i]);
    FreeParticleMemory(ParticleAcceleration[i]);
    FreeBaryonFieldMemory(AccelerationField[i]);
    delete [] RandomForcingField[i];
  }
 
  FreeParticleMemory(ParticleAcceleration[MAX_DIMENSION]);
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    FreeBaryonFieldMemory(BaryonField[i]);
    FreeBaryonFieldMemory(OldBaryonField[i]);
    FreeBaryonFieldMemory(InterpolatedField[i]);
  }
   delete [] YT_TemperatureField;

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++) {
    if(OldAccelerationField[i] != NULL ){
      FreeBaryonFieldMemory(OldAccelerationField[i]);
      OldAccelerationField[i] = NULL;
    }
  }
#endif
 
  DeleteFluxes(BoundaryFluxes);
  delete BoundaryFluxes;
 
  FreeParticleMemory(ParticleMass);
  FreeParticleMemory(ParticleNumber);
  FreeParticleMemory(ParticleType);
  FreeBaryonFieldMemory(PotentialField);
  FreeBaryonFieldMemory(GravitatingMassField);
  FreeBaryonFieldMemory(GravitatingMassFieldParticles);
  delete [] FlaggingField;
  FreeBaryonFieldMemory(MassFlaggingField);
  FreeBaryonFieldMemory(ParticleMassFlaggingField);
 
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    FreeParticleMemory(ParticleAttribute[i]);

  delete [] divB;
  for (int i=0; i<3; i++) {
    delete gradPhi[i];
  }


  DeleteStarList(Stars);

#ifdef TRANSFER
  delete PhotonPackages;
  if (FinishedPhotonPackages != NULL)
    delete FinishedPhotonPackages;
  if (PausedPhotonPackages != NULL)
    delete PausedPhotonPackages;
  delete [] SubgridMarker;
#endif

/* 
  if (debug && GridRank > 0) {
    printf("grid->destructor: deleting grid with dims = ");
    WriteListOfInts(stdout, GridRank, GridDimension);
  }
*/
 
}
