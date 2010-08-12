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
/  Linear Newton system setup function
/
/  written by: Daniel Reynolds
/  date:       August 2009
/  modified1:  
/
/  PURPOSE: Called by implicit solver to notify the Problem of 
/           updates to the current state (given in the vector u), 
/           so that the linear system matrix J(u) may be updated 
/           if necessary.  For the multi-frequency problem, we here 
/           compute only the local Jacobian components over the domain, 
/           and leave the actual matrix setup/solve for the lsolve 
/           routine, since we use a Schur-complement formulation for 
/           the linear system solution.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"



int MFProb::lsetup(EnzoVector *u)
{

  float stime = MPI_Wtime();

//   if (debug)  printf("Entering MFProb::lsetup routine\n");

  // check that the MFProb has been set up
  if (!prepared) 
    ENZO_FAIL("MFProb lsetup: MFProb not yet prepared!");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("MFProb lsetup: x0 vector dims do not match!");
  if (usz[1] != LocDims[1]) 
    ENZO_FAIL("MFProb lsetup: x1 vector dims do not match!");
  if (usz[2] != LocDims[2]) 
    ENZO_FAIL("MFProb lsetup: x2 vector dims do not match!");
  if (usz[3] != (4+Nchem)) 
    ENZO_FAIL("MFProb lsetup: nspecies do not match!");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("MFProb lsetup: x0 vector sizes do not match!");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("MFProb lsetup: x1 vector sizes do not match!");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("MFProb lsetup: x2 vector sizes do not match!");

  // clear Jacobian data arrays
  int i, j, k, ns, outidx;
  L_e->constant(0.0);
  L_HI->constant(0.0);
  L_HeI->constant(0.0);
  L_HeII->constant(0.0);
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E1_E1[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E1_HI[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E2_E2[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E2_HI[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E2_HeI[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E3_E3[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E3_HI[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E3_HeI[i] = 0.0;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  L_E3_HeII[i] = 0.0;
  

  // depending on 'Model' we either approximate or analytically compute the 
  // local Jacobian components
  //    analytically compute Jacobian  
  //    ['Models' with analytic Jacobian right now]
  if (approx_jac == 0) {

    // analytical Jacobians not implemented (due to radiation composition and 
    // numerical quadratures); issue warning and continue with approximate Jacobian 
    if (debug) {
      printf("MFProb lsetup warning: analytical Jacobians not implemented!\n");
      printf("   Continuing with approximate Jacobians.\n");
    }
  }

  //    approximate Jacobian
//  else {

    float stime1 = MPI_Wtime();

    // we first compute the local radiation Jacobian components analytically
    //   (automatically stores these in L_E*_* arrays)
    if (this->RadJac(u) == FAIL) 
      ENZO_FAIL("MFProb lsetup: RadJac failure!");

    float ftime1 = MPI_Wtime();
    timers[20] += ftime1 - stime1;

    // next, we get ready to perform finite-differences

    stime1 = MPI_Wtime();

    //      get typical values for input vectors
    float utypical[4+Nchem];
    utypical[0] = U0->rmsnorm_component(0);  // radiation energy
    utypical[1] = U0->rmsnorm_component(1);  // radiation energy
    utypical[2] = U0->rmsnorm_component(2);  // radiation energy
    for (i=4; i<Nchem+4; i++)                // chemistry 
      utypical[i] = U0->rmsnorm_component(i);
    //      override fluid energy correction typical value (since ec0=0)
    float dtmp1=0.0, dtmp2=0.0;
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      dtmp1 += eh[i]*eh[i];
    dtmp1 = sqrt(dtmp1/ArrDims[0]/ArrDims[1]/ArrDims[2])/eScale;
#ifdef USE_MPI
    if (layout[0]*layout[1]*layout[2] == 1)
      utypical[3] = dtmp1;
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg one = 1;
      MPI_Allreduce(&dtmp1, &dtmp2, one, DataType, MPI_SUM, MPI_COMM_WORLD);
      utypical[3] = dtmp2/NumberOfProcessors;  // estimate based on equidistribution
    }
#else
    utypical[3] = dtmp1;
#endif
    //      just in case any of the species are not currently used
    for (i=0; i<Nchem+4; i++) 
      utypical[i] = (utypical[i] == 0.0) ? 1.0 : utypical[i];

    ftime1 = MPI_Wtime();
    timers[21] += ftime1 - stime1;

//     printf("MFProb_lsetup: utypical =\n");
//     for (i=0; i<Nchem+4; i++)  printf("   %g\n",utypical[i]);


    //   set up temporary arrays
    EnzoVector *fval = tmp1;
    EnzoVector *utmp = tmp2;
    EnzoVector *ftmp = tmp3;

//     float dtmp;
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHI[i]*piHI[i];
//     printf("   ||piHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeI[i]*piHeI[i];
//     printf("   ||piHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeII[i]*piHeII[i];
//     printf("   ||piHeII|| = %g\n",sqrt(dtmp));
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHI[i]*GHI[i];
//     printf("   ||GHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeI[i]*GHeI[i];
//     printf("   ||GHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0; for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeII[i]*GHeII[i];
//     printf("   ||GHeII|| = %g\n",sqrt(dtmp));

    //   compute the local resid at the current state (fval)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(fval,u) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    }
    else {
      if (this->LocResid(fval,u) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    }

//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||f(u%"ISYM")|| = %g\n",i,fval->rmsnorm_component(i));


    //    determine floating-point roundoff
    float epsilon=1.0;
    while ((1.0 + epsilon*0.5) > 1.0)  epsilon*=0.5;
    epsilon = sqrt(epsilon);

    // copy u into utmp
    utmp->copy(u);

//     printf("\n Jacobians wrt E1\n");
  
    stime1 = MPI_Wtime();

    // next, we approximate the Jacobians wrt E1
    int iz, iy, ix, idx, ns2;
    float sigma, *uarray, *utmparray, *farray, *ftmparray, *Lblock;
    float typfac = 0.001;
    //    perturb E1
    uarray = u->GetData(iE1);
//     printf("   ||u(E1)|| = %g\n",u->rmsnorm_component(iE1));
    utmparray = utmp->GetData(iE1);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[0]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      utmparray[idx] += sigma;
    }
//     printf("   ||utmp(E1)|| = %g\n",utmp->rmsnorm_component(iE1));

    //    update the photo-ionization and photo-heating rates
    if (this->ComputeRadiationIntegrals(utmp) == FAIL) 
      ENZO_FAIL("MFProb lsetup: ComputeRadiationIntegrals failure!");
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHI[i]*piHI[i];
//     printf("   ||piHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeI[i]*piHeI[i];
//     printf("   ||piHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeII[i]*piHeII[i];
//     printf("   ||piHeII|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHI[i]*GHI[i];
//     printf("   ||GHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeI[i]*GHeI[i];
//     printf("   ||GHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeII[i]*GHeII[i];
//     printf("   ||GHeII|| = %g\n",sqrt(dtmp));
    //    compute the local resid due to this perturbation (ftmp)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    }
    else {
      if (this->LocResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    } 
    // restore utmp
    utmp->copy(u);

//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


    //    store the resulting Jacobian approximations
    //       Jec_E1
    farray = fval->GetData(iec);
    ftmparray = ftmp->GetData(iec);
    Lblock = L_e->GetData(iE1);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[0]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_e_E1|| = %g\n",L_e->rmsnorm_component(iE1));
    //       JHI_E1
    farray = fval->GetData(iHI);
    ftmparray = ftmp->GetData(iHI);
    Lblock = L_HI->GetData(iE1);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[0]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_HI_E1|| = %g\n",L_HI->rmsnorm_component(iE1));
    if (Nchem == 3) {
      //       JHeI_E1
      farray = fval->GetData(iHeI);
      ftmparray = ftmp->GetData(iHeI);
      Lblock = L_HeI->GetData(iE1);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[0]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//       printf("   ||L_HeI_E1|| = %g\n",L_HeI->rmsnorm_component(iE1));
      //       JHeII_E1
      farray = fval->GetData(iHeII);
      ftmparray = ftmp->GetData(iHeII);
      Lblock = L_HeII->GetData(iE1);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[0]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//       printf("   ||L_HeII_E1|| = %g\n",L_HeII->rmsnorm_component(iE1));
    }

    ftime1 = MPI_Wtime();
    timers[22] += ftime1 - stime1;


//     printf("\n Jacobians wrt E2\n");
  
    stime1 = MPI_Wtime();

    // next, we approximate the Jacobians wrt E2
    //    perturb E2
    uarray = u->GetData(iE2);
//     printf("   ||u(E2)|| = %g\n",u->rmsnorm_component(iE2));
    utmparray = utmp->GetData(iE2);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[1]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      utmparray[idx] += sigma;
    }
//     printf("   ||utmp(E2)|| = %g\n",utmp->rmsnorm_component(iE2));
    //    update the photo-ionization and photo-heating rates
    if (this->ComputeRadiationIntegrals(utmp) == FAIL) 
      ENZO_FAIL("MFProb lsetup: ComputeRadiationIntegrals failure!");
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHI[i]*piHI[i];
//     printf("   ||piHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeI[i]*piHeI[i];
//     printf("   ||piHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeII[i]*piHeII[i];
//     printf("   ||piHeII|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHI[i]*GHI[i];
//     printf("   ||GHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeI[i]*GHeI[i];
//     printf("   ||GHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeII[i]*GHeII[i];
//     printf("   ||GHeII|| = %g\n",sqrt(dtmp));
    //    compute the local resid due to this perturbation (ftmp)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    } 
    else {
      if (this->LocResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


    //    store the resulting Jacobian approximations
    //       Jec_E2
    farray = fval->GetData(iec);
    ftmparray = ftmp->GetData(iec);
    Lblock = L_e->GetData(iE2);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[1]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_e_E2|| = %g\n",L_e->rmsnorm_component(iE2));
    //       JHI_E2
    farray = fval->GetData(iHI);
    ftmparray = ftmp->GetData(iHI);
    Lblock = L_HI->GetData(iE2);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[1]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_HI_E2|| = %g\n",L_HI->rmsnorm_component(iE2));
    if (Nchem == 3) {
      //       JHeI_E2
      farray = fval->GetData(iHeI);
      ftmparray = ftmp->GetData(iHeI);
      Lblock = L_HeI->GetData(iE2);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[1]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeI_E2|| = %g\n",L_HeI->rmsnorm_component(iE2));
      //       JHeII_E2
      farray = fval->GetData(iHeII);
      ftmparray = ftmp->GetData(iHeII);
      Lblock = L_HeII->GetData(iE2);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[1]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeII_E2|| = %g\n",L_HeII->rmsnorm_component(iE2));
    }
    ftime1 = MPI_Wtime();
    timers[23] += ftime1 - stime1;


//     printf("\n Jacobians wrt E3\n");
  
    stime1 = MPI_Wtime();

    // next, we approximate the Jacobians wrt E3
    //    perturb E3
    uarray = u->GetData(iE3);
//     printf("   ||u(E3)|| = %g\n",u->rmsnorm_component(iE3));
    utmparray = utmp->GetData(iE3);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[2]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      utmparray[idx] += sigma;
    }
//     printf("   ||utmp(E3)|| = %g\n",utmp->rmsnorm_component(iE3));
    //    update the photo-ionization and photo-heating rates
    if (this->ComputeRadiationIntegrals(utmp) == FAIL) 
      ENZO_FAIL("MFProb lsetup: ComputeRadiationIntegrals failure!");
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHI[i]*piHI[i];
//     printf("   ||piHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeI[i]*piHeI[i];
//     printf("   ||piHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += piHeII[i]*piHeII[i];
//     printf("   ||piHeII|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHI[i]*GHI[i];
//     printf("   ||GHI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeI[i]*GHeI[i];
//     printf("   ||GHeI|| = %g\n",sqrt(dtmp));
//     dtmp=0.0;  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
//       dtmp += GHeII[i]*GHeII[i];
//     printf("   ||GHeII|| = %g\n",sqrt(dtmp));
    //    compute the local resid due to this perturbation (ftmp)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    } 
    else {
      if (this->LocResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


    //    store the resulting Jacobian approximations
    //       Jec_E3
    farray = fval->GetData(iec);
    ftmparray = ftmp->GetData(iec);
    Lblock = L_e->GetData(iE3);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[2]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_e_E3|| = %g\n",L_e->rmsnorm_component(iE3));
    //       JHI_E3
    farray = fval->GetData(iHI);
    ftmparray = ftmp->GetData(iHI);
    Lblock = L_HI->GetData(iE3);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[2]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_HI_E3|| = %g\n",L_HI->rmsnorm_component(iE3));
    if (Nchem == 3) {
      //       JHeI_E3
      farray = fval->GetData(iHeI);
      ftmparray = ftmp->GetData(iHeI);
      Lblock = L_HeI->GetData(iE3);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[2]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeI_E3|| = %g\n",L_HeI->rmsnorm_component(iE3));
      //       JHeII_E3
      farray = fval->GetData(iHeII);
      ftmparray = ftmp->GetData(iHeII);
      Lblock = L_HeII->GetData(iE3);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[2]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeII_E3|| = %g\n",L_HeII->rmsnorm_component(iE3));
    }

    ftime1 = MPI_Wtime();
    timers[24] += ftime1 - stime1;


    stime1 = MPI_Wtime();

    //    re-set the photo-ionization and photo-heating rates
    if (this->ComputeRadiationIntegrals(u) == FAIL) 
      ENZO_FAIL("MFProb lsetup: ComputeRadiationIntegrals failure!");

    ftime1 = MPI_Wtime();
    timers[25] += ftime1 - stime1;


//     printf("\n Jacobians wrt ec\n");

    stime1 = MPI_Wtime();
  
    // next, we approximate the Jacobians wrt ec
    //    perturb ec
    uarray = u->GetData(iec);
//     printf("   ||u(ec)|| = %g\n",u->rmsnorm_component(iec));
    utmparray = utmp->GetData(iec);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[3]);
      utmparray[idx] += sigma;
    }
    //    compute the local resid due to this perturbation (ftmp)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    } 
    else {
      if (this->LocResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));

   
    //    store the resulting Jacobian approximations
    //       Jec_ec
    farray = fval->GetData(iec);
    ftmparray = ftmp->GetData(iec);
    Lblock = L_e->GetData(iec);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[3]);
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_e_e|| = %g\n",L_e->rmsnorm_component(iec));
    //       JHI_ec
    farray = fval->GetData(iHI);
    ftmparray = ftmp->GetData(iHI);
    Lblock = L_HI->GetData(iec);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[3]);
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_HI_e|| = %g\n",L_HI->rmsnorm_component(iec));
    if (Nchem == 3) {
      //       JHeI_ec
      farray = fval->GetData(iHeI);
      ftmparray = ftmp->GetData(iHeI);
      Lblock = L_HeI->GetData(iec);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[3]);
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeI_e|| = %g\n",L_HeI->rmsnorm_component(iec));
      //       JHeII_ec
      farray = fval->GetData(iHeII);
      ftmparray = ftmp->GetData(iHeII);
      Lblock = L_HeII->GetData(iec);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[3]);
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeII_e|| = %g\n",L_HeII->rmsnorm_component(iec));
    }

    ftime1 = MPI_Wtime();
    timers[26] += ftime1 - stime1;

//     printf("\n Jacobians wrt HI\n");
  
    stime1 = MPI_Wtime();

    // next, we approximate the Jacobians wrt HI
    //    perturb HI
    uarray = u->GetData(iHI);
//     printf("   ||u(HI)|| = %g\n",u->rmsnorm_component(iHI));
    utmparray = utmp->GetData(iHI);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[4]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      utmparray[idx] += sigma;
    }
//     printf("   ||utmp(HI)|| = %g\n",utmp->rmsnorm_component(iHI));
    //    compute the local resid due to this perturbation (ftmp)
    if (AnalyticChem == 1) {
      if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
    } 
    else {
      if (this->LocResid(ftmp,utmp) == FAIL) 
	ENZO_FAIL("MFProb lsetup: LocResid failure!");
    } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


    //    store the resulting Jacobian approximations
    //       Jec_HI
    farray = fval->GetData(iec);
    ftmparray = ftmp->GetData(iec);
    Lblock = L_e->GetData(iHI);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[4]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_e_HI|| = %g\n",L_e->rmsnorm_component(iHI));
    //       JHI_HI
    farray = fval->GetData(iHI);
    ftmparray = ftmp->GetData(iHI);
    Lblock = L_HI->GetData(iHI);
    for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
      // sigma = epsilon*max(fabs(uarray[idx]),1.0);
      sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[4]);
      sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
      Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
    }
//     printf("   ||L_HI_HI|| = %g\n",L_HI->rmsnorm_component(iHI));
    if (Nchem == 3) {
      //       JHeI_HI
      farray = fval->GetData(iHeI);
      ftmparray = ftmp->GetData(iHeI);
      Lblock = L_HeI->GetData(iHI);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[4]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeI_HI|| = %g\n",L_HeI->rmsnorm_component(iHI));
      //       JHeII_HI
      farray = fval->GetData(iHeII);
      ftmparray = ftmp->GetData(iHeII);
      Lblock = L_HeII->GetData(iHI);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[4]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HeII_HI|| = %g\n",L_HeII->rmsnorm_component(iHI));
    }

    ftime1 = MPI_Wtime();
    timers[27] += ftime1 - stime1;

    if (Nchem == 3) {

//       printf("\n Jacobians wrt HeI\n");
  
      stime1 = MPI_Wtime();

      // next, we approximate the Jacobians wrt HeI
      //    perturb HeI
      uarray = u->GetData(iHeI);
//     printf("   ||u(HeI)|| = %g\n",u->rmsnorm_component(iHeI));
      utmparray = utmp->GetData(iHeI);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[5]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	utmparray[idx] += sigma;
      }
//     printf("   ||utmp(HeI)|| = %g\n",utmp->rmsnorm_component(iHeI));
      //    compute the local resid due to this perturbation (ftmp)
      if (AnalyticChem == 1) {
	if (this->AnalyticResid(ftmp,utmp) == FAIL)
	  ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
      } 
      else {
	if (this->LocResid(ftmp,utmp) == FAIL) 
	  ENZO_FAIL("MFProb lsetup: LocResid failure!");
      } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


      //    store the resulting Jacobian approximations
      //       Jec_HeI
      farray = fval->GetData(iec);
      ftmparray = ftmp->GetData(iec);
      Lblock = L_e->GetData(iHeI);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[5]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_e_HeI|| = %g\n",L_e->rmsnorm_component(iHeI));
      //       JHI_HeI
      farray = fval->GetData(iHI);
      ftmparray = ftmp->GetData(iHI);
      Lblock = L_HI->GetData(iHeI);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[5]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HI_HeI|| = %g\n",L_HI->rmsnorm_component(iHeI));
      if (Nchem == 3) {
	//       JHeI_HeI
	farray = fval->GetData(iHeI);
	ftmparray = ftmp->GetData(iHeI);
	Lblock = L_HeI->GetData(iHeI);
	for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	  // sigma = epsilon*max(fabs(uarray[idx]),1.0);
	  sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[5]);
	  sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	  Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
	}
//     printf("   ||L_HeI_HeI|| = %g\n",L_HeI->rmsnorm_component(iHeI));
	//       JHeII_HeI
	farray = fval->GetData(iHeII);
	ftmparray = ftmp->GetData(iHeII);
	Lblock = L_HeII->GetData(iHeI);
	for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	  // sigma = epsilon*max(fabs(uarray[idx]),1.0);
	  sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[5]);
	  sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	  Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
	}
//     printf("   ||L_HeII_HeI|| = %g\n",L_HeII->rmsnorm_component(iHeI));
      }

      ftime1 = MPI_Wtime();
      timers[28] += ftime1 - stime1;

//       printf("\n Jacobians wrt HeII\n");
  
      stime1 = MPI_Wtime();

      // next, we approximate the Jacobians wrt HeII
      //    perturb HeII
      uarray = u->GetData(iHeII);
//     printf("   ||u(HeII)|| = %g\n",u->rmsnorm_component(iHeII));
      utmparray = utmp->GetData(iHeII);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[6]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	utmparray[idx] += sigma;
      }
//     printf("   ||utmp(HeII)|| = %g\n",utmp->rmsnorm_component(iHeII));
      //    compute the local resid due to this perturbation (ftmp)
      if (AnalyticChem == 1) {
	if (this->AnalyticResid(ftmp,utmp) == FAIL) 
	  ENZO_FAIL("MFProb lsetup: AnalyticResid failure!");
      } 
      else {
	if (this->LocResid(ftmp,utmp) == FAIL) 
	  ENZO_FAIL("MFProb lsetup: LocResid failure!");
      } 
    // restore utmp
    utmp->copy(u);


//     for (i=0; i<Nchem+4; i++)  
//       printf("   ||ftmp(%"ISYM")|| = %g\n",i,ftmp->rmsnorm_component(i));


      //    store the resulting Jacobian approximations
      //       Jec_HeII
      farray = fval->GetData(iec);
      ftmparray = ftmp->GetData(iec);
      Lblock = L_e->GetData(iHeII);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[6]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_e_HeII|| = %g\n",L_e->rmsnorm_component(iHeII));
      //       JHI_HeII
      farray = fval->GetData(iHI);
      ftmparray = ftmp->GetData(iHI);
      Lblock = L_HI->GetData(iHeII);
      for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	// sigma = epsilon*max(fabs(uarray[idx]),1.0);
	sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[6]);
	sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
      }
//     printf("   ||L_HI_HeII|| = %g\n",L_HI->rmsnorm_component(iHeII));
      if (Nchem == 3) {
	//       JHeI_HeII
	farray = fval->GetData(iHeI);
	ftmparray = ftmp->GetData(iHeI);
	Lblock = L_HeI->GetData(iHeII);
	for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	  // sigma = epsilon*max(fabs(uarray[idx]),1.0);
	  sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[6]);
	  sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	  Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
	}
//     printf("   ||L_HeI_HeII|| = %g\n",L_HeI->rmsnorm_component(iHeII));
	//       JHeII_HeII
	farray = fval->GetData(iHeII);
	ftmparray = ftmp->GetData(iHeII);
	Lblock = L_HeII->GetData(iHeII);
	for (idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++) {
	  // sigma = epsilon*max(fabs(uarray[idx]),1.0);
	  sigma = epsilon*max(fabs(uarray[idx]),typfac*utypical[6]);
	  sigma *= (uarray[idx] > sigma) ? -1.0 : 1.0;
	  Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
	}
//     printf("   ||L_HeII_HeII|| = %g\n",L_HeII->rmsnorm_component(iHeII));
      }

      ftime1 = MPI_Wtime();
      timers[29] += ftime1 - stime1;
    }

//  }  // end approximate Jacobian

  // in case we are using the AnalyticChem solver with Model 4 (isothermal), 
  // we need to put something in the J_{ec,ec} block to avoid numerical 
  // problems (even though it is decoupled from the other equations).
  if ((AnalyticChem == 1) && (Model == 4))  L_e->addconst_component(iec,1.0);


  float ftime = MPI_Wtime();
  timers[13] += ftime - stime;

  // return success
  return SUCCESS;
}
#endif
