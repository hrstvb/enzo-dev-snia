/*****************************************************************************
 *                                                                           *
 * Copyright 2006 James Bordner
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Board of Trustees of the University of Illinois            *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
//======================================================================
//
// File:        performance.C
//
// Description: Performance-related code
//
//----------------------------------------------------------------------
//
// Namespaces:  jb
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

#include "performance.h"
#include "macros_and_parameters.h"

#ifdef USE_LCAPERF

void lcaperfInitialize (int max_level)
{

  // Initialize lcaperf

  lcaperf.initialize ("out.lcaperf");

  // Define lcaperf attributes

  lcaperf.new_attribute ("cycle", LCAP_INT);
  lcaperf.new_attribute ("level", LCAP_INT);

  // Define lcaperf counters

  lcaperf.new_counter ("count-zones",      user_counter_type_absolute);
  lcaperf.new_counter ("count-ghosts",     user_counter_type_absolute);
  lcaperf.new_counter ("count-grids",      user_counter_type_absolute);
  lcaperf.new_counter ("count-particles",  user_counter_type_absolute);
  lcaperf.new_counter ("time-sim",         user_counter_type_absolute);

  // Select which regions to print()

  lcaperf.new_region ("EvolveHierarchy");
  lcaperf.new_region ("EvolveLevel");
  lcaperf.new_region ("grid_SolveHydroEquations");
  lcaperf.new_region ("EvolvePhotons");
  lcaperf.new_region ("RebuildHierarchy");
  lcaperf.new_region ("SetBoundaryConditions");
  lcaperf.new_region ("ComovingExpansionTerms");
  lcaperf.new_region ("ComputePotentialFieldLevelZero");
  lcaperf.new_region ("ComputeRandomForcingNormalization");
  lcaperf.new_region ("ComputeRandomForcingNormalization");
  lcaperf.new_region ("ComputeRandomForcingNormalization");
  lcaperf.new_region ("CopyOverlappingMassField");
  lcaperf.new_region ("CreateFluxes");
  lcaperf.new_region ("DepositParticleMassField");
  lcaperf.new_region ("GetProjectedBoundaryFluxes");
  lcaperf.new_region ("grid_AddRandomForcing");
  lcaperf.new_region ("grid_ComputeAccelerationFieldExternal");
  lcaperf.new_region ("grid_CopyBaryonFieldToOldBaryonField");
  lcaperf.new_region ("grid_MultiSpeciesHandler");
  lcaperf.new_region ("grid_ParticleSplitter");
  lcaperf.new_region ("grid_SolveForPotential");
  lcaperf.new_region ("grid_StarParticleHandler");
  lcaperf.new_region ("InlineHaloFinder");
  lcaperf.new_region ("PrepareDensityField");
  lcaperf.new_region ("PrepareGravitatingMassField1");
  lcaperf.new_region ("PrepareGravitatingMassField2a");
  lcaperf.new_region ("PrepareGravitatingMassField2b");
  lcaperf.new_region ("ProjectSolutionToParentGrid");
  lcaperf.new_region ("RadiationFieldUpdate");
  lcaperf.new_region ("SetBC_Parent");
  lcaperf.new_region ("SetBC_Siblings");
  lcaperf.new_region ("SetLevelTimeStep");
  lcaperf.new_region ("SolveForPotential");
  lcaperf.new_region ("StarParticleAccretion");
  lcaperf.new_region ("StarParticleAddFeedback");
  lcaperf.new_region ("StarParticleDeath");
  lcaperf.new_region ("StarParticleFinalize");
  lcaperf.new_region ("StarParticleInitialize");
  lcaperf.new_region ("StarParticleSubtractAccretedMass");
  lcaperf.new_region ("star_FindFeedbackSphere");
  lcaperf.new_region ("star_FindFeedbackSphere2");
  lcaperf.new_region ("star_FindFeedbackSphere_Sum");
  lcaperf.new_region ("star_FindFeedbackSphere_Zero");
  lcaperf.new_region ("star_SphereContained");
  lcaperf.new_region ("star_UpdatePositionVelocity");
  lcaperf.new_region ("UpdateFromFinerGrids");
  lcaperf.new_region ("UpdateParticlePositions");
  lcaperf.new_region ("UpdateStarParticleCount");

  for (int level=0; level <= max_level; level++) {

    char lcaperf_counter_name[30];

    sprintf (lcaperf_counter_name,"count-zones-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-ghosts-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-grids-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-particles-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

  }

  lcaperf.begin();
}

//----------------------------------------------------------------------

void lcaperfFinalize ()
{
  lcaperf.end();
  lcaperf.finalize();
}

#endif /* USE_LCAPERF */


