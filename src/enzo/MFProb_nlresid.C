/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Frequency, Multi-species, Implicit Problem Class
/  Nonlinear residual function
/
/  written by: Daniel Reynolds
/  date:       August 2009
/  modified1:  
/
/  PURPOSE: Nonlinear residual function that defines the coupled,
/           implicit-time, radiation diffusion/chemistry/fluid energy 
/           system.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"



int MFProb::nlresid(EnzoVector *fu, EnzoVector *u)
{

  float stime = MPI_Wtime();
  
//   if (debug)
//     printf("Entering MFProb::nlresid routine\n");

  // check that the MFProb has been set up
  if (!prepared)
    ENZO_FAIL("MFProb nlresid: MFProb not yet prepared!");

//   // have u communicate neighbor information
//   if (u->exchange() == FAIL) 
//     ENZO_FAIL("MFProb nlresid: EnzoVector exchange failure!");
  // have u initiate communication of neighbor information
  if (u->exchange_start() == FAIL) 
    ENZO_FAIL("MFProb nlresid: EnzoVector exchange_start failure!");

  // initialize residual to zero
  if (fu->constant(0.0) == FAIL) 
    ENZO_FAIL("MFProb nlresid: EnzoVector constant failure!");
  
  // have u finish communication of neighbor information
  if (u->exchange_end() == FAIL) 
    ENZO_FAIL("MFProb nlresid: EnzoVector exchange_end failure!");

  // enforce boundary conditions on state u
  if (this->EnforceBoundary(u,0) == FAIL) 
    ENZO_FAIL("MFProb nlresid: EnforceBoundary failure!");

  //    update the photo-ionization and photo-heating rates
  if (this->ComputeRadiationIntegrals(u) == FAIL) 
    ENZO_FAIL("MFProb nlresid: ComputeRadiationIntegrals failure!");

  // compute the radiation residuals
  if (this->RadResid(fu, u) == FAIL) 
    ENZO_FAIL("MFProb nlresid: RadResid failure!");

  // compute the chemistry & gas energy residuals depending on AnalyticChem flag
  if (AnalyticChem == 1) {
    if (this->AnalyticResid(fu, u) == FAIL) 
      ENZO_FAIL("MFProb nlresid: AnalyticResid failure!");
  }
  else {
    if (this->LocResid(fu, u) == FAIL) 
      ENZO_FAIL("MFProb nlresid: LocResid failure!");
  }

  //   Enforce boundary conditions on fu (for Dirichlet faces)
  if (this->EnforceBoundary(fu,1) == FAIL) 
    ENZO_FAIL("MFProb nlresid: EnforceBoundary failure!");

  
//   for (int i=0; i<4+Nchem; i++) 
//     printf("   f(%"ISYM") = %g\n",i,fu->rmsnorm_component(i));


  float ftime = MPI_Wtime();
  timers[15] += ftime - stime;

  // return success
  return SUCCESS;

}
#endif
