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
/  Linear Newton system solution function
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/
/  PURPOSE: Solves the linear Newton system J(u)*s = b.  For the 
/           multi-frequency Problem without advection, the 
/           problem may be solved via the following Schur-complement 
/           formulation.
/
/           The Jacobian may be written in the operator form: 
/                   [ L_e_e  L_e_n      L_e_E     ]   
/               J = [ L_n_e  L_n_n      L_n_E     ] = [ M U ]
/                   [ L_E_e  L_E_n  (L_E_E + D_E) ]   [ L D ]
/           where for the fluid energy e, chemical species n, and 
/           radiation group E,
/               L_ee = local Jacobian of e wrt e
/               L_en = local Jacobian of e wrt n_j, j=1:Nchem
/               L_eE = local Jacobian of e wrt E_k, k=1:3
/               L_ne = local Jacobian of n_i wrt e, i=1:Nchem
/               L_nn = local Jacobian of n_i wrt n_j, i,j=1:Nchem
/               L_nE = local Jacobian of n_i wrt E_k, i=1:Nchem, k=1:3
/               L_Ee = local Jacobian of E_k wrt e, k=1:3
/               L_En = local Jacobian of E_k wrt n_j, j=1:Nchem, k=1:3
/               L_EE = local Jacobian of E_k wrt E_j, k=1:3, j=1:3
/               D_EE = diffusive Jacobians of E_k wrt E_k, k=1:3
/                 M = [ L_ee L_en ]
/                     [ L_ne L_nn ]
/                 U = [ L_eE L_nE ]^T
/                 L = [ L_Ee L_En ]
/                 D = (L_EE + D_EE).
/            The Schur complement formulation provides that
/                 Ji = [ I -Mi*U ][ Mi 0  ][   I   0 ]
/                      [ 0    I  ][ 0  Pi ][ -L*Mi I ]
/            where we use 'i' to denote the inverse, e.g. Ai = A^{-1}, 
/            and where the Schur complement is formed as P = D-L*Mi*U.  
/            Therefore, the solve J*s = b, where 
/            s = (s_e s_n s_E)^T = (s_m s_E)^T  (other vectors similarly 
/            condense e and n into m) may be broken down into the stages:
/                 (1) Solve for c_m:  M*c_m = b_m             (local)
/                 (2) Solve for Y_m:  M*Y_m = U               (local)
/                 (3) Update: b_E = b_E - L*c_m               (local)
/                 (4) Construct:  Y_E = L_EE - L*Y_m          (local)
/                 (5) Construct: P = D_EE + Y_E               (local)
/                 (6) Solve for s_E:  P*s_E = b_E             (nonlocal)
/                 (7) Set s_m = c_m - Y_m*s_E.                (local)
/            We note that all of the steps are completely local except 
/            for the step (6), which requires one vector-valued diffusive 
/            solve, for which we use the HYPRE library's PFMG solver.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"


int MFProb::lsolve(EnzoVector *s, EnzoVector *b, EnzoVector *u, float delta)
{
#ifdef USE_JBPERF
    JBPERF_START("MFProb_lsolve");
#endif
  float stime = MPI_Wtime();

//   if (debug)  printf("Entering MFProb::lsolve routine\n");

#ifdef USE_HYPRE

  // check that the MFProb has been set up
  if (!prepared) 
    ENZO_FAIL("MFProb lsolve: MFProb not yet prepared!");
  
  // have b communicate neighbor information and enforce BCs
//   if (b->exchange() == FAIL) 
//     ENZO_FAIL("MFProb lsolve: EnzoVector exchange failure!");
  if (this->EnforceBoundary(b,1) == FAIL) 
    ENZO_FAIL("MFProb lsolve: EnforceBoundary failure!");

  // check that b matches local vector size (including ghosts, etc)
  int ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  int vsz[4];
  b->size(&vsz[0], &vsz[1], &vsz[2], &vsz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (vsz[0] != LocDims[0]) 
    ENZO_FAIL("MFProb lsolve: x0 vector dims do not match!");
  if (vsz[1] != LocDims[1]) 
    ENZO_FAIL("MFProb lsolve: x1 vector dims do not match!");
  if (vsz[2] != LocDims[2]) 
    ENZO_FAIL("MFProb lsolve: x2 vector dims do not match!");
  if (vsz[3] != (4+Nchem)) 
    ENZO_FAIL("MFProb lsolve: nspecies do not match!");
  if ((vsz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("MFProb lsolve: x0 vector sizes do not match!");
  if ((vsz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("MFProb lsolve: x1 vector sizes do not match!");
  if ((vsz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("MFProb lsolve: x2 vector sizes do not match!");

//   if (debug)  printf("MFProb::lsolve -- performing steps 1-2\n");

  // clear temporary vectors for use in Schur complement correction
  tmp1->constant(0.0);
  tmp2->constant(0.0);
  tmp3->constant(0.0);

  float stime1 = MPI_Wtime();

  //////////////////////////////////////////////////////////////
  // steps (1) and (2): local solves 
  //          c_m = Mi*b_m
  //          y_m = Mi*U
  int ix, iy, iz, idx, size=Nchem+1, four=4;
  float *M = new float[size*size];
  float *bvec = new float[4*size];
  float *xvec = new float[4*size];
  //    set accessor pointers to relevant arrays
  float *YE1_E1, *YE2_E1, *YE3_E1, *Yec_E1, *Yn1_E1, *Yn2_E1, *Yn3_E1;
  float *YE1_E2, *YE2_E2, *YE3_E2, *Yec_E2, *Yn1_E2, *Yn2_E2, *Yn3_E2;
  float *YE1_E3, *YE2_E3, *YE3_E3, *Yec_E3, *Yn1_E3, *Yn2_E3, *Yn3_E3;
  float *L_ec_E1, *L_ec_E2, *L_ec_E3, *L_ec_ec, *L_ec_n1, *L_ec_n2, *L_ec_n3;
  float *L_n1_E1, *L_n1_E2, *L_n1_E3, *L_n1_ec, *L_n1_n1, *L_n1_n2, *L_n1_n3;
  float *L_n2_E1, *L_n2_E2, *L_n2_E3, *L_n2_ec, *L_n2_n1, *L_n2_n2, *L_n2_n3;
  float *L_n3_E1, *L_n3_E2, *L_n3_E3, *L_n3_ec, *L_n3_n1, *L_n3_n2, *L_n3_n3;
  float *b_E1, *b_E2, *b_E3, *b_ec, *b_n1, *b_n2, *b_n3;
  float *c_E1, *c_E2, *c_E3, *c_ec, *c_n1, *c_n2, *c_n3;
  YE1_E1 = tmp1->GetData(iE1);
  YE2_E1 = tmp1->GetData(iE2);
  YE3_E1 = tmp1->GetData(iE3);
  Yec_E1 = tmp1->GetData(iec);
  Yn1_E1 = tmp1->GetData(iHI);
  YE1_E2 = tmp2->GetData(iE1);
  YE2_E2 = tmp2->GetData(iE2);
  YE3_E2 = tmp2->GetData(iE3);
  Yec_E2 = tmp2->GetData(iec);
  Yn1_E2 = tmp2->GetData(iHI);
  YE1_E3 = tmp3->GetData(iE1);
  YE2_E3 = tmp3->GetData(iE2);
  YE3_E3 = tmp3->GetData(iE3);
  Yec_E3 = tmp3->GetData(iec);
  Yn1_E3 = tmp3->GetData(iHI);
  L_ec_E1 = L_e->GetData(iE1);
  L_ec_E2 = L_e->GetData(iE2);
  L_ec_E3 = L_e->GetData(iE3);
  L_ec_ec = L_e->GetData(iec);
  L_ec_n1 = L_e->GetData(iHI);
  L_n1_E1 = L_HI->GetData(iE1);
  L_n1_E2 = L_HI->GetData(iE2);
  L_n1_E3 = L_HI->GetData(iE3);
  L_n1_ec = L_HI->GetData(iec);
  L_n1_n1 = L_HI->GetData(iHI);
  b_E1 = b->GetData(iE1);
  b_E2 = b->GetData(iE2);
  b_E3 = b->GetData(iE3);
  b_ec = b->GetData(iec);
  b_n1 = b->GetData(iHI);
  if (Nchem == 3) {
    Yn2_E1 = tmp1->GetData(iHeI);
    Yn3_E1 = tmp1->GetData(iHeII);
    Yn2_E2 = tmp2->GetData(iHeI);
    Yn3_E2 = tmp2->GetData(iHeII);
    Yn2_E3 = tmp3->GetData(iHeI);
    Yn3_E3 = tmp3->GetData(iHeII);
    L_ec_n2 = L_e->GetData(iHeI);
    L_ec_n3 = L_e->GetData(iHeII);
    L_n1_n2 = L_HI->GetData(iHeI);
    L_n1_n3 = L_HI->GetData(iHeII);
    L_n2_E1 = L_HeI->GetData(iE1);
    L_n2_E2 = L_HeI->GetData(iE2);
    L_n2_E3 = L_HeI->GetData(iE3);
    L_n2_ec = L_HeI->GetData(iec);
    L_n2_n1 = L_HeI->GetData(iHI);
    L_n2_n2 = L_HeI->GetData(iHeI);
    L_n2_n3 = L_HeI->GetData(iHeII);
    L_n3_E1 = L_HeII->GetData(iE1);
    L_n3_E2 = L_HeII->GetData(iE2);
    L_n3_E3 = L_HeII->GetData(iE3);
    L_n3_ec = L_HeII->GetData(iec);
    L_n3_n1 = L_HeII->GetData(iHI);
    L_n3_n2 = L_HeII->GetData(iHeI);
    L_n3_n3 = L_HeII->GetData(iHeII);
    b_n2 = b->GetData(iHeI);
    b_n3 = b->GetData(iHeII);
  }

//   printf("MFProb_lsolve:\n");
//   printf("    ||L_ec_E1|| = %g\n",L_e->rmsnorm_component(iE1));
//   printf("    ||L_ec_E2|| = %g\n",L_e->rmsnorm_component(iE2));
//   printf("    ||L_ec_E3|| = %g\n",L_e->rmsnorm_component(iE3));
//   printf("    ||L_ec_ec|| = %g\n",L_e->rmsnorm_component(iec));
//   printf("    ||L_ec_n1|| = %g\n",L_e->rmsnorm_component(iHI));
//   printf("    ||L_n1_E1|| = %g\n",L_HI->rmsnorm_component(iE1));
//   printf("    ||L_n1_E2|| = %g\n",L_HI->rmsnorm_component(iE2));
//   printf("    ||L_n1_E3|| = %g\n",L_HI->rmsnorm_component(iE3));
//   printf("    ||L_n1_ec|| = %g\n",L_HI->rmsnorm_component(iec));
//   printf("    ||L_n1_n1|| = %g\n",L_HI->rmsnorm_component(iHI));
//   printf("    ||L_ec_n2|| = %g\n",L_e->rmsnorm_component(iHeI));
//   printf("    ||L_ec_n3|| = %g\n",L_e->rmsnorm_component(iHeII));
//   printf("    ||L_n1_n2|| = %g\n",L_HI->rmsnorm_component(iHeI));
//   printf("    ||L_n1_n3|| = %g\n",L_HI->rmsnorm_component(iHeII));
//   printf("    ||L_n2_E1|| = %g\n",L_HeI->rmsnorm_component(iE1));
//   printf("    ||L_n2_E2|| = %g\n",L_HeI->rmsnorm_component(iE2));
//   printf("    ||L_n2_E3|| = %g\n",L_HeI->rmsnorm_component(iE3));
//   printf("    ||L_n2_ec|| = %g\n",L_HeI->rmsnorm_component(iec));
//   printf("    ||L_n2_n1|| = %g\n",L_HeI->rmsnorm_component(iHI));
//   printf("    ||L_n2_n2|| = %g\n",L_HeI->rmsnorm_component(iHeI));
//   printf("    ||L_n2_n3|| = %g\n",L_HeI->rmsnorm_component(iHeII));
//   printf("    ||L_n3_E1|| = %g\n",L_HeII->rmsnorm_component(iE1));
//   printf("    ||L_n3_E2|| = %g\n",L_HeII->rmsnorm_component(iE2));
//   printf("    ||L_n3_E3|| = %g\n",L_HeII->rmsnorm_component(iE3));
//   printf("    ||L_n3_ec|| = %g\n",L_HeII->rmsnorm_component(iec));
//   printf("    ||L_n3_n1|| = %g\n",L_HeII->rmsnorm_component(iHI));
//   printf("    ||L_n3_n2|| = %g\n",L_HeII->rmsnorm_component(iHeI));
//   printf("    ||L_n3_n3|| = %g\n",L_HeII->rmsnorm_component(iHeII));
//   printf("    ||b_E1|| = %g\n",b->rmsnorm_component(iE1));
//   printf("    ||b_E2|| = %g\n",b->rmsnorm_component(iE2));
//   printf("    ||b_E3|| = %g\n",b->rmsnorm_component(iE3));
//   printf("    ||b_ec|| = %g\n",b->rmsnorm_component(iec));
//   printf("    ||b_n1|| = %g\n",b->rmsnorm_component(iHI));
//   printf("    ||b_n2|| = %g\n",b->rmsnorm_component(iHeI));
//   printf("    ||b_n3|| = %g\n",b->rmsnorm_component(iHeII));

  //          set shortcuts for c_m (store in b_m)
  c_E1 = b_E1;
  c_E2 = b_E2;
  c_E3 = b_E3;
  c_ec = b_ec;
  c_n1 = b_n1;
  c_n2 = b_n2;
  c_n3 = b_n3;

  //          perform computations depending on Nchem
  if (Nchem == 1) {
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;

	// set up local multiple RHS matrix system M*[c_m, Y_m] = [b_m, U]
	//      M matrix (column-major [Fortran] order)
	M[0] = L_ec_ec[idx];
	M[1] = L_n1_ec[idx];

	M[2] = L_ec_n1[idx];
	M[3] = L_n1_n1[idx];

	//      b_m (column-major [Fortran] order)
	bvec[0] = b_ec[idx];
	bvec[1] = b_n1[idx];

	//      U (column-major [Fortran] order)
	bvec[2] = L_ec_E1[idx];
	bvec[3] = L_n1_E1[idx];

	bvec[4] = L_ec_E2[idx];
	bvec[5] = L_n1_E2[idx];

	bvec[6] = L_ec_E3[idx];
	bvec[7] = L_n1_E3[idx];

	// solve dense local systems
	if (this->BlockSolve(M, xvec, bvec, &size, &four) != SUCCESS) 
	  ENZO_FAIL("MFProb lsolve: BlockSolve failure!");

	// extract solution components to appropriate locations
	//      c_m (column-major [Fortran] order)
	c_ec[idx] = xvec[0];
	c_n1[idx] = xvec[1];

	//      Y_m (column-major [Fortran] order)
	Yec_E1[idx] = xvec[2];
	Yn1_E1[idx] = xvec[3];

	Yec_E2[idx] = xvec[4];
	Yn1_E2[idx] = xvec[5];

	Yec_E3[idx] = xvec[6];
	Yn1_E3[idx] = xvec[7];
      } } }
  }
  else {
    for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;

	  // set up local multiple RHS matrix system M*[c_m, Y_m] = [b_m, U]
	  //      M matrix (column-major [Fortran] order)
	  M[0] = L_ec_ec[idx];
	  M[1] = L_n1_ec[idx];
	  M[2] = L_n2_ec[idx];
	  M[3] = L_n3_ec[idx];

	  M[4] = L_ec_n1[idx];
	  M[5] = L_n1_n1[idx];
	  M[6] = L_n2_n1[idx];
	  M[7] = L_n3_n1[idx];

	  M[8]  = L_ec_n2[idx];
	  M[9]  = L_n1_n2[idx];
	  M[10] = L_n2_n2[idx];
	  M[11] = L_n3_n2[idx];

	  M[12] = L_ec_n3[idx];
	  M[13] = L_n1_n3[idx];
	  M[14] = L_n2_n3[idx];
	  M[15] = L_n3_n3[idx];

	  //      b_m (column-major [Fortran] order)
	  bvec[0] = b_ec[idx];
	  bvec[1] = b_n1[idx];
	  bvec[2] = b_n2[idx];
	  bvec[3] = b_n3[idx];

	  //      U (column-major [Fortran] order)
	  bvec[4] = L_ec_E1[idx];
	  bvec[5] = L_n1_E1[idx];
	  bvec[6] = L_n2_E1[idx];
	  bvec[7] = L_n3_E1[idx];

	  bvec[8]  = L_ec_E2[idx];
	  bvec[9]  = L_n1_E2[idx];
	  bvec[10] = L_n2_E2[idx];
	  bvec[11] = L_n3_E2[idx];

	  bvec[12] = L_ec_E3[idx];
	  bvec[13] = L_n1_E3[idx];
	  bvec[14] = L_n2_E3[idx];
	  bvec[15] = L_n3_E3[idx];

	  // solve dense local systems
	  if (this->BlockSolve(M, xvec, bvec, &size, &four) != SUCCESS) 
	    ENZO_FAIL("MFProb lsolve: BlockSolve failure!");
	  
	  // extract solution components to appropriate locations
	  //      c_m (column-major [Fortran] order)
	  c_ec[idx] = xvec[0];
	  c_n1[idx] = xvec[1];
	  c_n2[idx] = xvec[2];
	  c_n3[idx] = xvec[3];

	  //      Y_m (column-major [Fortran] order)
	  Yec_E1[idx] = xvec[4];
 	  Yn1_E1[idx] = xvec[5];
	  Yn2_E1[idx] = xvec[6];
	  Yn3_E1[idx] = xvec[7];

	  Yec_E2[idx] = xvec[8];
	  Yn1_E2[idx] = xvec[9];
	  Yn2_E2[idx] = xvec[10];
	  Yn3_E2[idx] = xvec[11];

	  Yec_E3[idx] = xvec[12];
	  Yn1_E3[idx] = xvec[13];
	  Yn2_E3[idx] = xvec[14];
	  Yn3_E3[idx] = xvec[15];
	} } }
  }
  delete[] M;
  delete[] xvec;
  delete[] bvec;

  float ftime1 = MPI_Wtime();
  timers[16] += ftime1 - stime1;

//   if (debug)  printf("MFProb::lsolve -- performing steps 3-4\n");

  stime1 = MPI_Wtime();

  //////////////////////////////////////////////////////////////
  // steps (3) and (4): Compute rhs c_E and Construct local update for P
  //     (3) update:     c_E = b_E  - L*c_m
  //     (4) construct:  Y_E = L_EE - L*Y_m
  if (Nchem == 1) {
    for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	  
	  // The L matrix looks like
	  //    [ 0 L_E1_HI ]
	  //    [ 0 L_E2_HI ]
	  //    [ 0 L_E3_HI ]
	  //
	  // The c_E vector looks like
	  //    [ c_E1, c_E2, c_E3 ]^T
	  //
	  // The b_E vector looks like
	  //    [ b_E1, b_E2, b_E3 ]^T
	  //
	  // The c_m vector looks like
	  //    [ c_ec, c_n1 ]^T

	  // compute c_E = b_E - L*c_m
	  c_E1[idx] = b_E1[idx] - (L_E1_HI[idx]*c_n1[idx]);
	  c_E2[idx] = b_E2[idx] - (L_E2_HI[idx]*c_n1[idx]);
	  c_E3[idx] = b_E3[idx] - (L_E3_HI[idx]*c_n1[idx]);
	  
	  // L_EE looks like
	  //    [ L_E1_E1    0       0    ]
	  //    [    0    L_E2_E2    0    ]
	  //    [    0       0    L_E3_E3 ]
	  //
	  // The Y_m matrix looks like
	  //    [ Yec_E1  Yec_E2  Yec_E3 ]
	  //    [ Yn1_E1  Yn1_E2  Yn1_E3 ]
	  // 
	  // The Y_E matrix looks like
	  //    [ YE1_E1  YE1_E2  YE1_E3 ]
	  //    [ YE2_E1  YE2_E2  YE2_E3 ]
	  //    [ YE3_E1  YE3_E2  YE3_E3 ]

	  // update Y_E = L_EE - L*Y_m (this will be 3x3 dense)
	  YE1_E1[idx] = L_E1_E1[idx] - (L_E1_HI[idx]*Yn1_E1[idx]);
	  YE2_E1[idx] = -(L_E2_HI[idx]*Yn1_E1[idx]);
	  YE3_E1[idx] = -(L_E3_HI[idx]*Yn1_E1[idx]);
	  YE1_E2[idx] = -(L_E1_HI[idx]*Yn1_E2[idx]);
	  YE2_E2[idx] = L_E2_E2[idx] - (L_E2_HI[idx]*Yn1_E2[idx]);
	  YE3_E2[idx] = -(L_E3_HI[idx]*Yn1_E2[idx]);
	  YE1_E3[idx] = -(L_E1_HI[idx]*Yn1_E3[idx]);
	  YE2_E3[idx] = -(L_E2_HI[idx]*Yn1_E3[idx]);
	  YE3_E3[idx] = L_E3_E3[idx] - (L_E3_HI[idx]*Yn1_E3[idx]);
	  
	} } }
  }
  else {
    for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	  
	  // The L matrix looks like
	  //    [ 0 L_E1_HI    0         0     ]
	  //    [ 0 L_E2_HI L_E2_HeI     0     ]
	  //    [ 0 L_E3_HI L_E3_HeI L_E3_HeII ]
	  // 
	  // The c_E vector looks like
	  //    [ c_E1, c_E2, c_E3 ]^T
	  //
	  // The b_E vector looks like
	  //    [ b_E1, b_E2, b_E3 ]^T
	  //
	  // The c_m vector looks like
	  //    [ c_ec, c_n1, c_n2, c_n3 ]^T

	  // compute c_E = b_E - L*c_m
	  c_E1[idx] = b_E1[idx] - (L_E1_HI[idx]*c_n1[idx]);
	  c_E2[idx] = b_E2[idx] - (L_E2_HI[idx] *c_n1[idx]
				  +L_E2_HeI[idx]*c_n2[idx]);
	  c_E3[idx] = b_E3[idx] - (L_E3_HI[idx]  *c_n1[idx]
				  +L_E3_HeI[idx] *c_n2[idx]
				  +L_E3_HeII[idx]*c_n3[idx]);
	  
	  // L_EE looks like
	  //    [ L_E1_E1    0       0    ]
	  //    [    0    L_E2_E2    0    ]
	  //    [    0       0    L_E3_E3 ]
	  //
	  // The Y_m matrix looks like
	  //    [ Yec_E1  Yec_E2  Yec_E3 ]
	  //    [ Yn1_E1  Yn1_E2  Yn1_E3 ]
	  //    [ Yn2_E1  Yn2_E2  Yn2_E3 ]
	  //    [ Yn3_E1  Yn3_E2  Yn3_E3 ]
	  // 
	  // The Y_E matrix looks like
	  //    [ YE1_E1  YE1_E2  YE1_E3 ]
	  //    [ YE2_E1  YE2_E2  YE2_E3 ]
	  //    [ YE3_E1  YE3_E2  YE3_E3 ]
	  
	  // update Y_E = L_EE - L*Y_m (this will be 3x3 dense)
	  YE1_E1[idx] = L_E1_E1[idx] - (L_E1_HI[idx]*Yn1_E1[idx]);
	  YE2_E1[idx] = -(L_E2_HI[idx] *Yn1_E1[idx] 
			 +L_E2_HeI[idx]*Yn2_E1[idx]);
	  YE3_E1[idx] = -(L_E3_HI[idx]  *Yn1_E1[idx] 
			 +L_E3_HeI[idx] *Yn2_E1[idx] 
			 +L_E3_HeII[idx]*Yn3_E1[idx]);
	  YE1_E2[idx] = -(L_E1_HI[idx]*Yn1_E2[idx]);
	  YE2_E2[idx] = L_E2_E2[idx] - (L_E2_HI[idx] *Yn1_E2[idx]
				       +L_E2_HeI[idx]*Yn2_E2[idx]);
	  YE3_E2[idx] = -(L_E3_HI[idx]  *Yn1_E2[idx]
			 +L_E3_HeI[idx] *Yn2_E2[idx]
			 +L_E3_HeII[idx]*Yn3_E2[idx]);
	  YE1_E3[idx] = -(L_E1_HI[idx]*Yn1_E3[idx]);
	  YE2_E3[idx] = -(L_E2_HI[idx] *Yn1_E3[idx]
			 +L_E2_HeI[idx]*Yn2_E3[idx]);
	  YE3_E3[idx] = L_E3_E3[idx] - (L_E3_HI[idx]  *Yn1_E3[idx]
				       +L_E3_HeI[idx] *Yn2_E3[idx]
		 	  	       +L_E3_HeII[idx]*Yn3_E3[idx]);
	  
	} } }
  }
  ftime1 = MPI_Wtime();
  timers[17] += ftime1 - stime1;

//     printf("MFProb_lsolve:\n");
//     printf("     ||YE1_E1|| = %g\n",tmp1->rmsnorm_component(iE1));
//     printf("     ||YE2_E1|| = %g\n",tmp1->rmsnorm_component(iE2));
//     printf("     ||YE3_E1|| = %g\n",tmp1->rmsnorm_component(iE3));
//     printf("     ||Yec_E1|| = %g\n",tmp1->rmsnorm_component(iec));
//     printf("     ||Yn1_E1|| = %g\n",tmp1->rmsnorm_component(iHI));
//     printf("     ||YE1_E2|| = %g\n",tmp2->rmsnorm_component(iE1));
//     printf("     ||YE2_E2|| = %g\n",tmp2->rmsnorm_component(iE2));
//     printf("     ||YE3_E2|| = %g\n",tmp2->rmsnorm_component(iE3));
//     printf("     ||Yec_E2|| = %g\n",tmp2->rmsnorm_component(iec));
//     printf("     ||Yn1_E2|| = %g\n",tmp2->rmsnorm_component(iHI));
//     printf("     ||YE1_E3|| = %g\n",tmp3->rmsnorm_component(iE1));
//     printf("     ||YE2_E3|| = %g\n",tmp3->rmsnorm_component(iE2));
//     printf("     ||YE3_E3|| = %g\n",tmp3->rmsnorm_component(iE3));
//     printf("     ||Yec_E3|| = %g\n",tmp3->rmsnorm_component(iec));
//     printf("     ||Yn1_E3|| = %g\n",tmp3->rmsnorm_component(iHI));



//   if (debug)  printf("MFProb::lsolve -- performing step 5\n");

  stime1 = MPI_Wtime();

  //////////////////////////////////////////////////////////////
  // step (5): Construct (locally) the Schur complement matrix 
  //    P = D-L*Mi*U.  In the code, this corresponds to the 
  //    matrix  P = D_EE + Y_E

//   //       communicate yvec to spread local corrections to neighbors
//   if (tmp1->exchange() == FAIL) 
//     ENZO_FAIL("MFProb lsolve: EnzoVector exchane (tmp1) failure!");
//   if (tmp2->exchange() == FAIL) 
//     ENZO_FAIL("MFProb lsolve: EnzoVector exchane (tmp2) failure!");
//   if (tmp3->exchange() == FAIL) 
//     ENZO_FAIL("MFProb lsolve: EnzoVector exchane (tmp3) failure!");

  //       set matrix values over grid
  float *u_E = u->GetData(0);
  float *u0_E = U0->GetData(0);
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  if (this->SetupSystem(P11tmpvec, P12tmpvec, P13tmpvec, P21tmpvec, P22tmpvec, 
			P23tmpvec, P31tmpvec, P32tmpvec, P33tmpvec, c_E1, 
			c_E2, c_E3, u, YE1_E1, YE1_E2, YE1_E3, YE2_E1, 
			YE2_E2, YE2_E3, YE3_E1, YE3_E2, YE3_E3) != SUCCESS) 
    ENZO_FAIL("MFProb lsolve: SetupSystem failure!");
  Eint32 zed=0;
  Eint32 one=1;
  Eint32 Sentries[7] = {0, 1, 2, 3, 4, 5, 6};
  Eint32 Nentries[1] = {stSize};
  Eint32 Lentries[1] = {stSize+1};
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,0,stSize,Sentries,P11tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,0,1,Nentries,P12tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,0,1,Lentries,P13tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,1,stSize,Sentries,P22tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,1,1,Nentries,P21tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,1,1,Lentries,P23tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,2,stSize,Sentries,P33tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,2,1,Nentries,P31tmpvec); 
  HYPRE_SStructMatrixSetBoxValues(P,0,ilower,iupper,2,1,Lentries,P32tmpvec); 


  //       assemble matrix
  HYPRE_SStructMatrixAssemble(P);

  ftime1 = MPI_Wtime();
  timers[18] += ftime1 - stime1;

//   if (debug)  printf("MFProb::lsolve -- performing step 6\n");

  //////////////////////////////////////////////////////////////
  // step (6): Solve the (nonlocal) Schur complement system,
  // i.e. solve for s_E:  P*s_E = c_E

  stime1 = MPI_Wtime();

  //       re-scale delta to relative residual and not actual
  delta /= b->rmsnorm();
  delta = min(delta, 1.0e-6);

  //       insert rhs vector into HYPRE vector rhsvec
  int Zbl, Ybl;
  int xBuff = ghXl-SolvOff[0];
  int yBuff = (ghYl-SolvOff[1])-SolvIndices[1][0];
  int zBuff = (ghZl-SolvOff[2])-SolvIndices[2][0];
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = c_E1[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(rhsvec, 0, ilower, iupper, 0, HYPREbuff);
    }
  }
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = c_E2[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(rhsvec, 0, ilower, iupper, 1, HYPREbuff);
    }
  }
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = c_E3[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(rhsvec, 0, ilower, iupper, 2, HYPREbuff);
    }
  }


  //       insert sol vector into HYPRE vector solvec
  ilower[1] = SolvIndices[1][0];  ilower[2] = SolvIndices[2][0];
  iupper[1] = SolvIndices[1][1];  iupper[2] = SolvIndices[2][1];
  size = (SolvIndices[0][1]-SolvIndices[0][0]+1)
        *(SolvIndices[1][1]-SolvIndices[1][0]+1)
        *(SolvIndices[2][1]-SolvIndices[2][0]+1);
  for (ix=0; ix<size; ix++)  P11tmpvec[ix] = 0.0;
  HYPRE_SStructVectorSetBoxValues(solvec, 0, ilower, iupper, 0, P11tmpvec);
  HYPRE_SStructVectorSetBoxValues(solvec, 0, ilower, iupper, 1, P11tmpvec);
  HYPRE_SStructVectorSetBoxValues(solvec, 0, ilower, iupper, 2, P11tmpvec);

  //       assemble vectors
//   if (debug)  printf("lsolve: calling HYPRE_SStructVectorAssemble\n");
  HYPRE_SStructVectorAssemble(solvec);
  HYPRE_SStructVectorAssemble(rhsvec);



//   if (debug)  printf("Writing out matrix to file P.mat\n");
//   HYPRE_SStructMatrixPrint("P.mat",P,0);

//   if (debug)  printf("Writing out rhs to file b.vec\n");
//   HYPRE_SStructVectorPrint("b.vec",rhsvec,0);


  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_SStructSolver solver;
  HYPRE_SStructSolver preconditioner;
  HYPRE_SStructPCGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &preconditioner);

  //          set preconditioner options
  HYPRE_SStructSysPFMGSetMaxIter(preconditioner, sol_maxit/4);
  HYPRE_SStructSysPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_SStructSysPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(preconditioner, sol_npost);
  HYPRE_SStructSysPFMGSetPrintLevel(preconditioner, sol_printl);
  HYPRE_SStructSysPFMGSetLogging(preconditioner, sol_log);

  //          set solver options
  if (rank > 1) {
    HYPRE_SStructPCGSetMaxIter(solver, sol_maxit);
    HYPRE_SStructPCGSetPrecond(solver, 
		     (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSysPFMGSolve,  
		     (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSysPFMGSetup, 
		      preconditioner);
  }
  else {    // ignore pfmg preconditioner for 1D tests (bug); increase CG its
    HYPRE_SStructPCGSetMaxIter(solver, sol_maxit*500);
  }
  if (delta != 0.0) {
//     HYPRE_SStructPCGSetAbsoluteTol(solver, Eflt64(delta));
    HYPRE_SStructPCGSetTol(solver, Eflt64(delta));
  }
  HYPRE_SStructPCGSetup(solver, P, rhsvec, solvec);

  //       solve the linear system
  HYPRE_SStructPCGSolve(solver, P, rhsvec, solvec);

  //       extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;
  Eint32 Sits=0;
  Eint32 Pits=0;
  HYPRE_SStructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_SStructPCGGetNumIterations(solver, &Sits);
  HYPRE_SStructSysPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug)
    printf("   HYPRE resid = %g (tol = %g), PCG = %i, PFMG = %i\n",
	   finalresid,delta,Sits,Pits);

//   if (debug)  printf("Writing out solution vector to file x.vec\n");
//   HYPRE_SStructVectorPrint("x.vec",solvec,0);

  //       extract values from solution vector 
  HYPRE_SStructVectorGather(solvec);
  //          s_E1
  float *s_E1 = s->GetData(iE1);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_SStructVectorGetBoxValues(solvec, 0, ilower, iupper, 0, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E1[Zbl+Ybl+xBuff+ix] = HYPREbuff[ix];
    }
  }
  //          s_E2
  float *s_E2 = s->GetData(iE2);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_SStructVectorGetBoxValues(solvec, 0, ilower, iupper, 1, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E2[Zbl+Ybl+xBuff+ix] = HYPREbuff[ix];
    }
  }
  //          s_E3
  float *s_E3 = s->GetData(iE3);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_SStructVectorGetBoxValues(solvec, 0, ilower, iupper, 2, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E3[Zbl+Ybl+xBuff+ix] = HYPREbuff[ix];
    }
  }

  //       destroy HYPRE matrix, vector and solver structures
  HYPRE_SStructPCGDestroy(solver);
  HYPRE_SStructSysPFMGDestroy(preconditioner);
  ftime1 = MPI_Wtime();
  timers[19] += ftime1 - stime1;

//   if (debug)  printf("MFProb::lsolve -- performing step 7\n");

  //////////////////////////////////////////////////////////////
  // step (7) Set s_m = c_m - Y_m*s_E
  float *s_ec = s->GetData(iec);
  float *s_n1 = s->GetData(iHI);
  if (Nchem == 1) {
    for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	  // The Y_m matrix looks like
	  //    [ Yec_E1    Yec_E2    Yec_E3   ]
	  //    [ Yn1_E1    Yn1_E2    Yn1_E3   ]
	  // 
	  // The s_E vector looks like
	  //    [ s_E1, s_E2, s_E3 ]^T
	  //
	  // The c_m vector looks like
	  //    [ c_ec, c_n1 ]^T
	  //
	  // The s_m vector looks like
	  //    [ s_ec, s_n1 ]^T
	  s_ec[idx] = c_ec[idx] - (Yec_E1[idx]*s_E1[idx]
				  +Yec_E2[idx]*s_E2[idx]
				  +Yec_E3[idx]*s_E3[idx]);
	  s_n1[idx] = c_n1[idx] - (Yn1_E1[idx]*s_E1[idx]
				  +Yn1_E2[idx]*s_E2[idx]
				  +Yn1_E3[idx]*s_E3[idx]);
	} } }
  } else{
    float *s_n2 = s->GetData(iHeI);
    float *s_n3 = s->GetData(iHeII);
    for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	  // The Y_m matrix looks like
	  //    [ Yec_E1  Yec_E2  Yec_E3 ]
	  //    [ Yn1_E1  Yn1_E2  Yn1_E3 ]
	  //    [ Yn2_E1  Yn2_E2  Yn2_E3 ]
	  //    [ Yn3_E1  Yn3_E2  Yn3_E3 ]
	  // 
	  // The s_E vector looks like
	  //    [ s_E1, s_E2, s_E3 ]^T
	  //
	  // The c_m vector looks like
	  //    [ c_ec, c_n1, c_n2, c_n3 ]^T
	  //
	  // The s_m vector looks like
	  //    [ s_ec, s_n1, s_n2, s_n3 ]^T
	  s_ec[idx] = c_ec[idx] - (Yec_E1[idx]*s_E1[idx]
				  +Yec_E2[idx]*s_E2[idx]
				  +Yec_E3[idx]*s_E3[idx]);
	  s_n1[idx] = c_n1[idx] - (Yn1_E1[idx]*s_E1[idx]
				  +Yn1_E2[idx]*s_E2[idx]
				  +Yn1_E3[idx]*s_E3[idx]);
	  s_n2[idx] = c_n2[idx] - (Yn2_E1[idx]*s_E1[idx]
			  	  +Yn2_E2[idx]*s_E2[idx]
				  +Yn2_E3[idx]*s_E3[idx]);
	  s_n3[idx] = c_n3[idx] - (Yn3_E1[idx]*s_E1[idx]
				  +Yn3_E2[idx]*s_E2[idx]
				  +Yn3_E3[idx]*s_E3[idx]);
	} } }
  }

//   if (debug)  printf("MFProb::lsolve -- finished!\n");

//   /////// TEMPORARY ///////
//     ENZO_FAIL("MFProb lsolve: diagnostics");
//   /////////////////////////


  float ftime = MPI_Wtime();
  timers[14] += ftime - stime;

#else

  ENZO_FAIL("MFProb lsolve Error: HYPRE must be enabled to use this module!");

#endif  // USE_HYPRE

  // return success
  return SUCCESS;
}
#endif
