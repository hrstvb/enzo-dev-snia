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
/             else -> analytical solution to local reaction network
/
************************************************************************/
#ifdef TRANSFER
#include "MFSplit.h"

 
 

int MFSplit::InitialGuess(EnzoVector *uvec)
{
  float stime = MPI_Wtime();

  // set initial guess depending on user choice
  switch (initial_guess) {


  case 0:  // previous time step

    uvec->copy(U0);
    break;

  case 1:
  case 2:
  case 3:
  case 4:
  case 5:  // analytical solution to local reaction network

    //   Note: extsrc is available thanks to the last call to nlresid
    uvec->copy(U0);
    if (this->AnalyticInitGuess(uvec) == FAIL)
      ENZO_FAIL("MFSplit InitialGuess: AnalyticInitGuess failure");
    break;

  default:  // illegal choice

    ENZO_VFAIL("MFSplit InitialGuess: illegal initial_guess choice = %"ISYM"\n",
	       initial_guess)

  }

  float ftime = MPI_Wtime();
  timers[12] += ftime - stime;

  return SUCCESS;

}

#endif
