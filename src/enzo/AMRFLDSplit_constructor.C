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
/  Split Implicit Problem Class, Constructor routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"


AMRFLDSplit::AMRFLDSplit()
{

//   if (debug)  printf("\nEntering AMRFLDSplit::constructor routine\n");
  int dim, face;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif

  // initialize total RT time to zero
  RTtime = 0.0;
  AMRSolTime = 0.0;

  // initialize solver parameters
  sol_tolerance = -1.0;
  sol_maxit = -1;
  sol_rlxtype = -1;
  sol_npre = -1;
  sol_npost = -1;
  sol_printl = -1;
  sol_log = -1;
  totIters = -1;

  // initialize problem grid information to -1/NULL
  rank = -1;
  for (dim=0; dim<3; dim++) {
    LocDims[dim] = 1;
    for (face=0; face<2; face++) {
      OnBdry[dim][face] = false;
      BdryType[dim][face] = -1;
      BdryVals[dim][face] = NULL;
    }
  }
  
  // initialize time-stepping related data to -1/NULL
  maxdt = 1.0e20;
  mindt = 0.0;
  initdt = 1.0e20;
  maxsubcycles = 1.0;
  dtfac = 1.0e20;
  dtnorm = 0.0;
  tnew = -1.0;
  told = -1.0;
  dt = -1.0;
  dtrad = -1.0;
  theta = -1.0;
  
  // initialize problem defining data 
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;
  aUnits = 1.0;
  ErScale = 1.0;
  ErUnits = 1.0;
  ErUnits0 = 1.0;
  NiUnits = 1.0;
  NiUnits0 = 1.0;
  DenUnits = 1.0;
  DenUnits0 = 1.0;
  LenUnits = 1.0;
  LenUnits0 = 1.0;
  TimeUnits = 1.0;
  VelUnits = 1.0;
  Nchem = 0;
  ESpectrum = -1;
  intSigE = 0.0;
  intSigESigHI = 0.0;
  intSigESigHeI = 0.0;
  intSigESigHeII = 0.0;
  intSigESigHInu = 0.0;
  intSigESigHeInu = 0.0;
  intSigESigHeIInu = 0.0;

}
#endif
