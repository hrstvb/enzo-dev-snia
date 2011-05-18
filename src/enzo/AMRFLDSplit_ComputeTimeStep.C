/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Computes the radiation time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
 
 
float AMRFLDSplit::ComputeTimeStep(Eflt64 Echange)
{

  // Set time step depending on how it has been set up by the user:
  //    If dtfac is set, compute maximum time step as estimate 
  //    allowing dtfac relative change.  This relative change is
  //    calculated elsewhere
  float dt_est = huge_number;    // max time step (normalized units)
  if (dtfac != huge_number) {

    // compute time step estimate (physical units)
    if (Echange <= 0.0) {
      dt_est = huge_number;
    } else {
      dt_est = dt*dtfac/Echange;
    }
    dt_est = min(dt_est, huge_number);

  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;
}

#endif
