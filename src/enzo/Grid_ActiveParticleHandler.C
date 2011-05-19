/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF ACTIVE PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  April, 2009 by JHW to have multiple types of star 
/              particles
/  modified2:  May, 2011 by MJT to be in support of active particles
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
int grid::ActiveParticleHandler(HierarchyEntry* SubgridPointer, int level,
                                float dtLevelAbove)
{

  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);

  /* initialize */

  LCAPERF_START("grid_ActiveParticleHandler");

    /* This is where everything used to be! */

  /* First we identify the data dependencies */

  struct ActiveParticleFormationDataFlags flags = flags_default;

  int i;
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleTypeToEvaluate->describe_data_flags(flags);
  }

  struct ActiveParticleFormationData supplemental_data;
  supplemental_data.level = level;
  std::vector<ActiveParticleType> temp_particles;
  supplemental_data.particles = &temp_particles;

  ActiveParticleType::ConstructData(this, flags, supplemental_data);

  std::vector<ActiveParticleType> NewParticles;

  int NumberOfNewParticles = 0;
  /* Now we iterate */
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    NumberOfNewParticles += ActiveParticleTypeToEvaluate->formation_function(
                                this, supplemental_data);
    NewParticles.insert(NewParticles.end(),
                        supplemental_data.particles->begin(),
                        supplemental_data.particles->end());
    supplemental_data.particles->clear();
  }

  ActiveParticleType::DestroyData(this, supplemental_data);

  /* Now we copy the particles from NewParticles into a statically allocated
   * array */

  int OldNumberOfActiveParticles = this->NumberOfActiveParticles;
  ActiveParticleType *OldActiveParticles = this->ActiveParticles;

  this->NumberOfActiveParticles += NumberOfNewParticles;
  this->ActiveParticles = new ActiveParticleType[this->NumberOfActiveParticles];
  for (i = 0; i < OldNumberOfActiveParticles; i++) {
    this->ActiveParticles[i] = OldActiveParticles[i];
  }
  for (i = 0; i < NumberOfNewParticles; i++) {
    this->ActiveParticles[this->NumberOfActiveParticles - i] =
        NewParticles.back();
    NewParticles.pop_back();
  }
  NewParticles.clear();
  //if (debug) printf("StarParticle: end\n");


  LCAPERF_STOP("grid_ActiveParticleHandler");
  return SUCCESS;
}
