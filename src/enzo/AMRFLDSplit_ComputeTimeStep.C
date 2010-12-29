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
/  PURPOSE: Computes the rad-hydro time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

 
 

float AMRFLDSplit::ComputeTimeStep(EnzoVector *uold, EnzoVector *unew)
{
  // get local mesh description
  int Nx, Ny, Nz, Nvar, ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  unew->size(&Nx, &Ny, &Nz, &Nvar, &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (Nx != LocDims[0]) 
    ENZO_FAIL("ComputeTimeStep error: x0 vector dims do not match");
  if (Ny != LocDims[1]) 
    ENZO_FAIL("ComputeTimeStep error: x1 vector dims do not match");
  if (Nz != LocDims[2]) 
    ENZO_FAIL("ComputeTimeStep error: x2 vector dims do not match");
  if (Nvar != 1) 
    ENZO_FAIL("ComputeTimeStep error: nspecies dims do not match");
  if ((Nx+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("ComputeTimeStep error: x0 vector sizes do not match");
  if ((Ny+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("ComputeTimeStep error: x1 vector sizes do not match");
  if ((Nz+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("ComputeTimeStep error: x2 vector sizes do not match");


  // Set time step depending on how it has been set up by the user:
  //    If dtfac is set, compute maximum time step as estimate 
  //    allowing dtfac relative change.  This relative change is
  //    estimated as follows:
  //       dt_new = dt_old / relerr_fac
  //    where relerr_fac gives the ratio between an estimated 
  //    local truncation error and the desired relative change:
  //       relerr_fac = || (unew - uold) / w ||_p
  //    with the scaling vector w given by
  //       w = dtfac*[sqrt(|unew*uold|) + atol]
  //    and where we have the following parameters:
  //       p - norm choice (input), 0->max norm, otherwise the p-norm
  //           **all p-norms here divide by the total number of cells**
  //       dtfac - desired relative change per step (input)
  //       atol - 0.1 (assumes units are all normalized)
  float dt_est = huge_number;    // max time step (normalized units)
  if (dtfac != huge_number) {

    // initialize variables
    float diff, w, tmp, atol;
    int i, j, k, l;
    int x0len = Nx + ghXl + ghXr;
    int x1len = Ny + ghYl + ghYr;
    float loc_est;

    // perform local estimates for the radiation energy relative change
    loc_est = 0.0;
    float *Eold = uold->GetData(0);
    float *Enew = unew->GetData(0);
    atol = 0.1; // assumes values are normalized
    if (dtnorm > 0.0) {
      for (k=ghZl; k<Nz+ghZl; k++) 
	for (j=ghYl; j<Ny+ghYl; j++)
	  for (i=ghXl; i<Nx+ghXl; i++) {
	    w = dtfac*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				*Eold[(k*x1len + j)*x0len + i])) 
		       + atol);
	    diff = Enew[(k*x1len + j)*x0len + i] 
	         - Eold[(k*x1len + j)*x0len + i];
	    tmp = fabs(diff/w);
	    loc_est += POW(tmp,dtnorm);
	  }
    }
    else {
      for (k=ghZl; k<Nz+ghZl; k++) 
	for (j=ghYl; j<Ny+ghYl; j++)
	  for (i=ghXl; i<Nx+ghXl; i++) {
	    w = dtfac*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				*Eold[(k*x1len + j)*x0len + i])) 
		       + atol);
	    diff = Enew[(k*x1len + j)*x0len + i] 
	         - Eold[(k*x1len + j)*x0len + i];
	    tmp = fabs(diff/w);
	    loc_est = (loc_est > tmp) ? loc_est : tmp;
	  }
    }

    // communicate to obtain overall sum/max
    float glob_est;
    int Nglobal = GlobDims[0]*GlobDims[1]*GlobDims[2];
#ifdef USE_MPI
    if (Nglobal == Nx*Ny*Nz)  glob_est = loc_est;
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg vars = 1;
      if (dtnorm > 0.0) 
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_SUM,MPI_COMM_WORLD);
      else
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_MAX,MPI_COMM_WORLD);
    }
#else
    glob_est = loc_est;
#endif

    // compute overall norms
    if (dtnorm > 0.0) {
      glob_est /= Nglobal;
      glob_est = POW(glob_est,1.0/dtnorm);
    }

    // compute time step estimate (physical units)
    dt_est = (glob_est == 0.0) ? huge_number : dt/glob_est;
    dt_est = min(dt_est, huge_number);
  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;

}

#endif
