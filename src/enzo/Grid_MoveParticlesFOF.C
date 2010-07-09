/***********************************************************************
/
/  GRID CLASS (MOVE ALL PARTICLES FROM SPECIFIED GRID TO FOF STRUCTURE)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#include "FOF_allvars.h"
 
int grid::MoveParticlesFOF(int level, FOF_particle_data* &P, 
			   int &Index, FOFData AllVars, float VelocityUnits, 
			   double MassUnits, int CopyDirection)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int i, j, dim, Index0, idx;
  double VelocityConv, MassConv;
  double VelConvInv, MassConvInv;
  double BoxSizeInv;

  VelocityConv = VelocityUnits / AllVars.UnitVelocity_in_cm_per_s;
  MassConv = MassUnits / pow(8.0, level) / AllVars.UnitMass_in_g;

  VelConvInv = 1.0 / VelocityConv;
  MassConvInv = 1.0 / MassConv;
  BoxSizeInv = 1.0 / AllVars.BoxSize;

  if (CopyDirection == COPY_OUT) {

    Index0 = Index;

#pragma omp parallel private(dim,j)
    {

    for (dim = 0; dim < GridRank; dim++) {
#pragma omp for nowait schedule(static) private(idx)
      for (i = 0; i < NumberOfParticles; i++) {
	idx = Index0+i;
	P[idx].Pos[dim] = AllVars.BoxSize * ParticlePosition[dim][i];
	P[idx].Vel[dim] = ParticleVelocity[dim][i] * VelocityConv;
      }
    }

    for (j = 0; j < NumberOfParticleAttributes; j++) {
#pragma omp for nowait schedule(static) 
      for (i = 0; i < NumberOfParticles; i++)
    	P[Index0+i].Attr[j] = ParticleAttribute[j][i];
    }

#pragma omp for nowait schedule(static) private(idx)
    for (i = 0; i < NumberOfParticles; i++) {
      idx = Index0+i;
      P[idx].slab = (int) (ParticlePosition[0][i] * NumberOfProcessors);

      P[idx].Mass = ParticleMass[i] * MassConv;

      P[idx].PartID = ParticleNumber[i];
      P[idx].Type = ParticleType[i];
      //P[idx].level = level;
      //P[idx].GridID = GridNum;

      P[idx].Energy = 0.0;
      P[idx].Rho = 0.0;

    }
    } // END omp parallel

    Index = Index0 + NumberOfParticles;
    this->DeleteParticles();

  } // ENDIF (COPY_OUT)

  else {

    /* When we create the smoothed DM fields, we copy particles across
       periodic boundaries.  Don't copy back these duplicated
       particles. */

    int npart = AllVars.Nlocal;//slab[MyProcessorNumber];
    for (i = 0; i < AllVars.Nlocal; i++)
      if (P[i].PartID < 0 || P[i].slab != MyProcessorNumber) npart--;

    /* Only move particles in this slab -- exclude the shadows. */

    this->AllocateNewParticles(npart);

#pragma omp parallel private(dim,j)
    {

    for (dim = 0; dim < GridRank; dim++)
#pragma omp for nowait schedule(static)
      for (i = 0; i < AllVars.Nlocal; i++)
	if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0) {
	  ParticlePosition[dim][i] = P[i].Pos[dim] * BoxSizeInv;
	  ParticleVelocity[dim][i] = P[i].Vel[dim] * VelConvInv;
	}

    for (j = 0; j < NumberOfParticleAttributes; j++)
#pragma omp for nowait schedule(static)
      for (i = 0; i < AllVars.Nlocal; i++)
	if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0)
	  ParticleAttribute[j][i] = P[i].Attr[j];

#pragma omp for nowait schedule(static)
    for (i = 0; i < AllVars.Nlocal; i++) 
      if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0) {
	ParticleMass[i] = P[i].Mass * MassConvInv;
	ParticleType[i] = P[i].Type;
	ParticleNumber[i] = P[i].PartID;
      }
    } // END omp parallel

    NumberOfParticles = npart;

  } // ENDELSE (COPY_IN)

  return SUCCESS;

}
