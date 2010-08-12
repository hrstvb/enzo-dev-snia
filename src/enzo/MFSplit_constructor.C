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
/  Multi-Frequency, Multi-species, Split Problem Class
/  Constructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "MFSplit.h"


MFSplit::MFSplit()
{

//   if (debug)  printf("\nEntering MFSplit::constructor routine\n");
  int dim, face;

  // initialize timer to 0
  RTtime = 0.0;
  for (int i=0; i<30; i++)  timers[i] = 0.0;

  // initialize free-streaming solver and array to NULL
  FSSolve = NULL;
  Efree = NULL;

  // initialize HYPRE values to -1/NULL
  stSize = -1;
#ifdef USE_HYPRE
  grid = NULL;
  stencil = NULL;
#endif
  sol_tolerance = -1.0;
  sol_maxit = -1;
  sol_rlxtype = -1;
  sol_npre = -1;
  sol_npost = -1;
  sol_printl = -1;
  sol_log = -1;
  totIters = -1;
  matentries = NULL;
  rhsentries = NULL;
  HYPREbuff = NULL;
  for (dim=0; dim<3; dim++) {
    for (face=0; face<2; face++)
      SolvIndices[dim][face] = 0;
    SolvOff[dim] = 0;
  }

  // initialize Newton solver values to -1/NULL
  initial_guess = -1;

  // initialize problem grid information to -1/NULL
  rank = -1;
  for (dim=0; dim<3; dim++) {
    layout[dim] = 1;    // initialize for single-processor run
    location[dim] = 0;  // initialize for single-processor run
    LocDims[dim] = 1;
    ArrDims[dim] = 1;
    dx[dim] = 1.0;
    for (face=0; face<2; face++) {
      OnBdry[dim][face] = false;
      NBors[dim][face] = MPI_PROC_NULL;
      GhDims[dim][face] = 0;
      BdryType[dim][face] = -1;
      EdgeVals[dim][face] = -1.0;
      EBdryVals[dim][face] = NULL;
    }
  }
  
  // initialize time-stepping related data to -1/NULL
  initdt = 1.0e20;
  maxdt = 1.0e20;
  mindt = 0.0;
  dtfac[0] = 1.0e20;
  dtfac[1] = 1.0e20;
  dtfac[2] = 1.0e20;
  dtnorm = 0.0;
  tnew = -1.0;
  told = -1.0;
  dt = -1.0;
  dtchem = -1.0;
  theta = -1.0;
  LimType = -1;
  sol = NULL;
  U0 = NULL;
  extsrc = NULL;

  // initialize problem defining data 
  iE1 = -1;
  iE2 = -1;
  iE3 = -1;
  iec = -1;
  iHI = -1;
  iHeI = -1;
  iHeII = -1;
  Nchem = 0;
  Model = -1;
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;
  aUnits = 1.0;
  rScale = 1.0;
  eScale = 1.0;
  nScale = 1.0;
  fsUnits = 1.0;
  rUnits = 1.0;
  rUnits0 = 1.0;
  E1Units = 1.0;
  E1Units0 = 1.0;
  E2Units = 1.0;
  E2Units0 = 1.0;
  E3Units = 1.0;
  E3Units0 = 1.0;
  eUnits = 1.0;
  nUnits = 1.0;
  nUnits0 = 1.0;
  DenUnits = 1.0;
  DenUnits0 = 1.0;
  LenUnits = 1.0;
  LenUnits0 = 1.0;
  TimeUnits = 1.0;
  VelUnits = 1.0;
  HFrac = 1.0;  // default is pure hydrogen (no helium)
  ESpectrum = -1;
  hnu0_HI = 0.0;
  hnu0_HeI = 0.0;
  hnu0_HeII = 0.0;
  piHI = NULL;
  piHeI = NULL;
  piHeII = NULL;
  GHI = NULL;
  GHeI = NULL;
  GHeII = NULL;
  chibar = 0.0;
  chinuint = 0.0;

  // initialize HYPRE structures to NULL
#ifdef USE_HYPRE
  P      = NULL;
  rhsvec = NULL;
  solvec = NULL;
#endif

  // initialize access to Enzo arrays to NULL
  vx = NULL;
  vy = NULL;
  vz = NULL;
  rho = NULL;
  eh = NULL;

  // initialize extra storage arrays to NULL
  eCorr = NULL;

}
#endif
