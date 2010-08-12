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
/  Initial Guess Computation Routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Computes an initial guess to the time-evolved solution, 
/           according to the method specified by the initial_guess
/           input parameter:
/                0 -> use previous time step
/                1 -> full fwd Euler (all rhs)
/                2 -> analytical solution to local reaction network
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"

 
 

int MFProb::InitialGuess(EnzoVector *uvec)
{

  float stime = MPI_Wtime();

  // set initial guess depending on user choice
  switch (initial_guess) {


  case 0:  // previous time step

    uvec->copy(U0);
    break;


  case 1:  // full forward Euler step

    if (this->nlresid(uvec, U0) == FAIL) 
      ENZO_FAIL("MFProb InitialGuess: nlresid failure!\n");

    // uguess = U0 - f(U0)  [f(U0) held in 'uvec']
    uvec->linearsum(1.0, U0, -1.0, uvec);
    break;


  case 2:  // analytical solution to local reaction network

    //   Note: extsrc is available thanks to the last call to nlresid
    uvec->copy(U0);
    if (this->AnalyticInitGuess(uvec) == FAIL) 
      ENZO_FAIL("MFProb InitialGuess: AnalyticInitGuess failure!\n");
    break;

  default:  // illegal choice

    ENZO_VFAIL("MFProb InitialGuess: illegal initial_guess choice = %"ISYM"\n",
	       initial_guess)

  }

  float ftime = MPI_Wtime();
  timers[12] += ftime - stime;

  return SUCCESS;
}

#endif
