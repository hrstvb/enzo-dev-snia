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
/  Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Computes the module time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called before MFProb::Return, so it 
/           sees whatever units that are internal to the solver module.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"

 
 

float MFProb::ComputeTimeStep(EnzoVector *unew)
{

  float stime = MPI_Wtime();

  // get local mesh description
  int Nx, Ny, Nz, Nvar, ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  unew->size(&Nx, &Ny, &Nz, &Nvar, &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (Nx != LocDims[0]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x0 vector dims do not match!\n");
  if (Ny != LocDims[1]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x1 vector dims do not match!\n");
  if (Nz != LocDims[2]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x2 vector dims do not match!\n");
  if (Nvar != (4+Nchem)) 
    ENZO_FAIL("MFProb ComputeTimeStep: nspecies dims do not match!\n");
  if ((Nx+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x0 vector sizes do not match!\n");
  if ((Ny+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x1 vector sizes do not match!\n");
  if ((Nz+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("MFProb ComputeTimeStep: x2 vector sizes do not match!\n");


  // Set time step depending on how it has been set up by the user:
  //    If dtfac is set for any species, compute maximum time step 
  //    as estimate allowing dtfac relative change.  This relative 
  //    change is estimated as follows:
  //       dt_new = dt_old / relerr_fac
  //    where relerr_fac gives the ratio between an estimated 
  //    local truncation error and the desired relative change:
  //       relerr_fac = || (unew - U0) / w ||_p
  //    with the scaling vector w given by
  //       w = dtfac*[sqrt(|unew*U0|) + atol]
  //    and where we have the following parameters:
  //       p - norm choice (input), 0->max norm, otherwise the p-norm
  //           **all p-norms here divide by the total number of cells**
  //       dtfac - desired relative change per step (input)
  //       atol - 1e-3 (assumes units are all normalized)
  //    For the gas energy correction, this is different since we do 
  //    not have U0: 
  //       relerr_fac = || unew / w ||_p
  //       w = dtfac*(unew + atol).
  float dt_est = maxdt;    // max time step (normalized units)
  float test = dtfac[0];
  int i, j, k, l, idx;
  for (i=0; i<4+Nchem; i++)  test = min(dtfac[i],test);
  if (test != huge_number) {

    // initialize variables
    float diff, w, tmp, atol;
    int x0len = Nx + ghXl + ghXr;
    int x1len = Ny + ghYl + ghYr;
    float loc_est[5];

    // perform local estimates for the radiation energy freq 1 relative change
    loc_est[0] = 0.0;
    if (dtfac[0] != huge_number) {
      float *Eold = U0->GetData(iE1);
      float *Enew = unew->GetData(iE1);
      atol = 0.01; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[0] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[0] = (loc_est[0] > tmp) ? loc_est[0] : tmp;
	    }
      }
    }

    // perform local estimates for the radiation energy freq 2 relative change
    loc_est[1] = 0.0;
    if (dtfac[0] != huge_number) {
      float *Eold = U0->GetData(iE2);
      float *Enew = unew->GetData(iE2);
      atol = 0.01; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[1] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[1] = (loc_est[1] > tmp) ? loc_est[1] : tmp;
	    }
      }
    }

    // perform local estimates for the radiation energy freq 3 relative change
    loc_est[2] = 0.0;
    if (dtfac[0] != huge_number) {
      float *Eold = U0->GetData(iE3);
      float *Enew = unew->GetData(iE3);
      atol = 0.01; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[2] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[2] = (loc_est[2] > tmp) ? loc_est[2] : tmp;
	    }
      }
    }

    // perform estimates for the gas energy
    loc_est[3] = 0.0;
    if (dtfac[1] != huge_number) {
      float *ec = unew->GetData(iec);
      atol = 0.01; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[1]*(fabs(ec[(k*x1len + j)*x0len + i]
			       + eh[(k*x1len + j)*x0len + i]/eScale) 
			    + atol);
	      diff = ec[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[3] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[1]*(fabs(ec[(k*x1len + j)*x0len + i]
			       + eh[(k*x1len + j)*x0len + i]/eScale) 
			    + atol);
	      diff = ec[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[3] = (loc_est[3] > tmp) ? loc_est[3] : tmp;
	    }
      }
    }

    // perform estimates for the chemistry
    float *HIold, *HInew, *HeIold, *HeInew, *HeIIold, *HeIInew;
    float ne_new, ne_old, nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII;
    float mp = 1.67262171e-24;
    loc_est[4] = 0.0;
    if (Nchem == 1) {
      if (dtfac[2] != huge_number) {
	HIold = U0->GetData(iHI);
	HInew = unew->GetData(iHI);
	atol = 0.01; // assumes values are normalized
	if (dtnorm > 0.0) {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		idx = (k*x1len + j)*x0len + i;
		nH = HFrac*rho[idx]*DenUnits/mp;
		nHI = HInew[idx]*nUnits;
		nHII = max(nH-nHI, 0.0);
		ne_new = nHII/nUnits;
		nH = HFrac*rho[idx]*DenUnits0/mp;
		nHI = HIold[idx]*nUnits0;
		nHII = max(nH-nHI, 0.0);
		ne_old = nHII/nUnits0;
		w = dtfac[2]*(sqrt(fabs(ne_new*ne_old))	+ atol);
		diff = ne_new - ne_old;
		tmp = fabs(diff/w);
		loc_est[4] += POW(tmp,dtnorm);
	      }
	}
	else {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		idx = (k*x1len + j)*x0len + i;
		nH = HFrac*rho[idx]*DenUnits/mp;
		nHI = HInew[idx]*nUnits;
		nHII = max(nH-nHI, 0.0);
		ne_new = nHII/nUnits;
		nH = HFrac*rho[idx]*DenUnits0/mp;
		nHI = HIold[idx]*nUnits0;
		nHII = max(nH-nHI, 0.0);
		ne_old = nHII/nUnits0;
		w = dtfac[2]*(sqrt(fabs(ne_new*ne_old))	+ atol);
		diff = ne_new - ne_old;
		tmp = fabs(diff/w);
		loc_est[4] = (loc_est[4] > tmp) ? loc_est[4] : tmp;
	      }
	}
      }
    } else {
	HIold = U0->GetData(iHI);
	HeIold = U0->GetData(iHeI);
	HeIIold = U0->GetData(iHeII);
	HInew = unew->GetData(iHI);
	HeInew = unew->GetData(iHeI);
	HeIInew = unew->GetData(iHeII);
	atol = 0.01; // assumes values are normalized
	if (dtnorm > 0.0) {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		idx = (k*x1len + j)*x0len + i;
		nH = HFrac*rho[idx]*DenUnits/mp;
		nHI = HInew[idx]*nUnits;
		nHII = max(nH-nHI, 0.0);
		nHe = (1.0-HFrac)*rho[idx]*DenUnits/mp;
		nHeI = HeInew[idx]*nUnits;
		nHeII = HeIInew[idx]*nUnits;
		nHeIII = max(nHe - nHeI - nHeII, 0.0);
		ne_new = (nHII + 0.25*nHeII + 0.5*nHeIII)/nUnits;
		nH = HFrac*rho[idx]*DenUnits0/mp;
		nHI = HIold[idx]*nUnits0;
		nHII = max(nH-nHI, 0.0);
		nHe = (1.0-HFrac)*rho[idx]*DenUnits0/mp;
		nHeI = HeIold[idx]*nUnits0;
		nHeII = HeIIold[idx]*nUnits0;
		nHeIII = max(nHe - nHeI - nHeII, 0.0);
		ne_old = (nHII + 0.25*nHeII + 0.5*nHeIII)/nUnits0;
		w = dtfac[2]*(sqrt(fabs(ne_new*ne_old))	+ atol);
		diff = ne_new - ne_old;
		tmp = fabs(diff/w);
		loc_est[4] += POW(tmp,dtnorm);
	      }
	}
	else {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		idx = (k*x1len + j)*x0len + i;
		nH = HFrac*rho[idx]*DenUnits/mp;
		nHI = HInew[idx]*nUnits;
		nHII = max(nH-nHI, 0.0);
		nHe = (1.0-HFrac)*rho[idx]*DenUnits/mp;
		nHeI = HeInew[idx]*nUnits;
		nHeII = HeIInew[idx]*nUnits;
		nHeIII = max(nHe - nHeI - nHeII, 0.0);
		ne_new = (nHII + 0.25*nHeII + 0.5*nHeIII)/nUnits;
		nH = HFrac*rho[idx]*DenUnits0/mp;
		nHI = HIold[idx]*nUnits0;
		nHII = max(nH-nHI, 0.0);
		nHe = (1.0-HFrac)*rho[idx]*DenUnits0/mp;
		nHeI = HeIold[idx]*nUnits0;
		nHeII = HeIIold[idx]*nUnits0;
		nHeIII = max(nHe - nHeI - nHeII, 0.0);
		ne_old = (nHII + 0.25*nHeII + 0.5*nHeIII)/nUnits0;
		w = dtfac[2]*(sqrt(fabs(ne_new*ne_old))	+ atol);
		diff = ne_new - ne_old;
		tmp = fabs(diff/w);
		loc_est[4] = (loc_est[4] > tmp) ? loc_est[4] : tmp;
	      }
	}
    }

    // communicate to obtain overall sum/max
    float glob_est[5];
    int Nglobal = GlobDims[0]*GlobDims[1]*GlobDims[2];
#ifdef USE_MPI
    if (Nglobal == Nx*Ny*Nz) 
      for (l=0; l<5; l++)  glob_est[l] = loc_est[l];
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg vars = 5;
      if (dtnorm > 0.0) 
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_SUM,MPI_COMM_WORLD);
      else
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_MAX,MPI_COMM_WORLD);
    }
#else
    for (l=0; l<5; l++)  glob_est[l] = loc_est[l];
#endif

    // compute overall norms
    if (dtnorm > 0.0) 
      for (l=0; l<5; l++) {
	glob_est[l] /= Nglobal;
	glob_est[l] = POW(glob_est[l],1.0/dtnorm);
      }

    // compute variable-specific time step estimates (physical units)
    float dt_est_var[5];
    for (l=0; l<5; l++) {
      dt_est_var[l] = (glob_est[l] == 0.0) ? huge_number : dt/glob_est[l];
      dt_est_var[l] = min(dt_est_var[l], huge_number);
    }

    // set estimated time step as minimum of component time steps
    dt_est = maxdt*TimeUnits;    // max time step estimate (physical units)
    for (l=0; l<5; l++) {
      dt_est = min(dt_est, dt_est_var[l]);
    }

    // limit maximum growth per step
    dt_est = min(dt_est, 1.1*dt);    // time step growth (physical units)

    // rescale dt estimates to normalized values
    dt_est /= TimeUnits;
    for (l=0; l<5; l++)  dt_est_var[l] /= TimeUnits;

    // account for min/max time step size (according to user)
    dt_est = max(dt_est, mindt);
    dt_est = min(dt_est, maxdt);

    if (debug) {
      printf("MFProb_ComputeTimestep: (E1,E2,E3,e,ne), dt = (");
      printf("%7.1e %7.1e %7.1e %7.1e %7.1e", dt_est_var[0], 
	     dt_est_var[1], dt_est_var[2], dt_est_var[3], dt_est_var[4]);
      printf(")\nMFProb_ComputeTimeStep: dt_est = %g\n",dt_est);
    }
  }

  float ftime = MPI_Wtime();
  timers[11] += ftime - stime;

  return dt_est;
}

#endif
