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
/  Linear Newton system solution function
/
/  written by: Daniel Reynolds
/  date:       December 2009
/
/  PURPOSE: Solves the linear radiation systems for the 
/           multi-frequency problem without advection, using the 
/           EnzoVector sol for both the initial guess, and for the final 
/           solution.  For this problem, we use the HYPRE library's 
/           system PFMG solver.
/
************************************************************************/
#ifdef TRANSFER
#include "MFSplit.h"


int MFSplit::LinearSolve(EnzoVector *sol, EnzoVector *src, float thisdt)
{
  float stime = MPI_Wtime();

//   if (debug)  printf("Entering MFSplit::LinearSolve routine\n");

#ifdef USE_HYPRE
  /////////// E1 equation ///////////
  //       set matrix values over grid
  float rhsnorm;
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, 
			1, sol, src, thisdt) != SUCCESS)
    ENZO_FAIL("MFSplit LinearSolve: SetupSystem failure (E1)");
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 

  //       assemble matrix
  HYPRE_StructMatrixAssemble(P);

  //       insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);

  //       set the linear solver tolerance (rescale to relative residual and not actual)
  Eflt64 delta = min(sol_tolerance/rhsnorm, 1.0e-6);

  //       insert zero initial guess into HYPRE vector solvec
  int ix, iy, iz, size;
  size = (SolvIndices[0][1]-SolvIndices[0][0]+1)
        *(SolvIndices[1][1]-SolvIndices[1][0]+1)
        *(SolvIndices[2][1]-SolvIndices[2][0]+1);
  for (ix=0; ix<size; ix++)  matentries[ix] = 0.0;
  HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, matentries);

  //       assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

//   if (debug)  printf("Writing out matrix to file P1.mat\n");
//   HYPRE_StructMatrixPrint("P1.mat",P,0);

//   if (debug)  printf("Writing out rhs to file b1.vec\n");
//   HYPRE_StructVectorPrint("b1.vec",rhsvec,0);

  //       for periodic dims, only coarsen until grid no longer divisible by 2
  Eint32 max_levels, level=-1;
  int Ndir;
  if (BdryType[0][0] == 0) {
    level = 0;
    Ndir = GlobDims[0];
    while ( Ndir%2 == 0 ) {
      level++;
      Ndir /= 2;
    }
  }
  max_levels = level;
  if (rank > 1) {
    if (BdryType[1][0] == 0) {
      level = 0;
      Ndir = GlobDims[0];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  if (rank > 2) {
    if (BdryType[2][0] == 0) {
      level = 0;
      Ndir = GlobDims[0];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  max_levels = min(level,max_levels);

  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_StructSolver solver1;
  HYPRE_StructSolver preconditioner1;
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver1);
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner1);

  //          set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner1, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner1, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner1, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner1, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner1, sol_npost);

  //          set solver options
  HYPRE_StructPCGSetPrintLevel(solver1, sol_printl);
  HYPRE_StructPCGSetLogging(solver1, sol_log);
  HYPRE_StructPCGSetRelChange(solver1, 1);
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver1, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver1, 
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
		      preconditioner1);
  }
  else {    // ignore pfmg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver1, sol_maxit*500);
  }
  if (delta != 0.0)  HYPRE_StructPCGSetTol(solver1, Eflt64(delta));
  HYPRE_StructPCGSetup(solver1, P, rhsvec, solvec);

  //       solve the linear system
  if (debug)  printf("  Monochromatic solve for E1:\n");
  HYPRE_StructPCGSolve(solver1, P, rhsvec, solvec);

//   if (debug)  printf("Writing out solution vector to file x1.vec\n");
//   HYPRE_StructVectorPrint("x1.vec",solvec,0);

  //       extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;
  Eint32 Sits=0;
  Eint32 Pits=0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver1, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver1, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner1, &Pits);
  totIters += Sits;
  if (debug)
    printf("    lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
	   finalresid*rhsnorm, sol_tolerance, rhsnorm, Sits, Pits);

  //       extract values from solution vector 
  int Zbl, Ybl, xBuff, yBuff, zBuff;
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  float *s_E = sol->GetData(iE1);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
    }
  }
  float ENorm1 = sol->rmsnorm_component(iE1);
  float EMax1  = sol->infnorm_component(iE1);
  if (debug) 
    printf("    E1 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   ENorm1,ENorm1*E1Units,EMax1,EMax1*E1Units);

  //       destroy HYPRE solver and preconditioner structures
  HYPRE_StructPCGDestroy(solver1);
  HYPRE_StructPFMGDestroy(preconditioner1);

  //       check for solver failure
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
    // if the final residual is too large, or is nan, quit
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {
      fprintf(stderr,"   Error: could not achieve prescribed tolerance!\n");

      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P1.mat\n");
      HYPRE_StructMatrixPrint("P1.mat",P,0);
      if (debug)  printf("Writing out rhs to file b1.vec\n");
      HYPRE_StructVectorPrint("b1.vec",rhsvec,0);
      if (debug)  printf("Writing out current solution to file x1.vec\n");
      HYPRE_StructVectorPrint("x1.vec",solvec,0);

      // dump module parameters to disk
      this->Dump(sol);

      fprintf(stderr,"  =====================================================================\n");
      ENZO_FAIL("MFSplit LinearSolve: HYPRE solver error (E1)");
    }
  }
  if (debug)  
    printf("  ---------------------------------------------------------------------\n");



  /////////// E2 equation ///////////
  //       set matrix values over grid
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, 
			2, sol, src, thisdt) != SUCCESS)
    ENZO_FAIL("MFSplit LinearSolve: SetupSystem failure (E2)");
  ilower[0] = SolvIndices[0][0];  iupper[0] = SolvIndices[0][1];
  ilower[1] = SolvIndices[1][0];  iupper[1] = SolvIndices[1][1];
  ilower[2] = SolvIndices[2][0];  iupper[2] = SolvIndices[2][1];
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 

  //       assemble matrix
  HYPRE_StructMatrixAssemble(P);

  //       insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);

  //       set the linear solver tolerance (rescale to relative residual and not actual)
  delta = min(sol_tolerance/rhsnorm, 1.0e-6);

  //       insert zero initial guess into HYPRE vector solvec
  for (ix=0; ix<size; ix++)  matentries[ix] = 0.0;
  HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, matentries);

  //       assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

//   if (debug)  printf("Writing out matrix to file P2.mat\n");
//   HYPRE_StructMatrixPrint("P2.mat",P,0);

//   if (debug)  printf("Writing out rhs to file b2.vec\n");
//   HYPRE_StructVectorPrint("b2.vec",rhsvec,0);

  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_StructSolver solver2;
  HYPRE_StructSolver preconditioner2;
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver2);
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner2);

  //          set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner2, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner2, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner2, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner2, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner2, sol_npost);

  //          set solver options
  HYPRE_StructPCGSetPrintLevel(solver2, sol_printl);
  HYPRE_StructPCGSetLogging(solver2, sol_log);
  HYPRE_StructPCGSetRelChange(solver2, 1);
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver2, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver2, 
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
		      preconditioner2);
  }
  else {    // ignore pfmg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver2, sol_maxit*500);
  }
  if (delta != 0.0)  HYPRE_StructPCGSetTol(solver2, Eflt64(delta));
  HYPRE_StructPCGSetup(solver2, P, rhsvec, solvec);

  //       solve the linear system
  if (debug)  printf("  Monochromatic solve for E2:\n");
  HYPRE_StructPCGSolve(solver2, P, rhsvec, solvec);

//   if (debug)  printf("Writing out solution vector to file x2.vec\n");
//   HYPRE_StructVectorPrint("x2.vec",solvec,0);

  //       extract solver & preconditioner statistics
  finalresid=1.0;
  Sits=0;
  Pits=0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver2, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver2, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner2, &Pits);
  totIters += Sits;
  if (debug)
    printf("    lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
	   finalresid*rhsnorm, sol_tolerance, rhsnorm, Sits, Pits);

  //       extract values from solution vector 
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  s_E = sol->GetData(iE2);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
    }
  }
  ENorm1 = sol->rmsnorm_component(iE2);
  EMax1  = sol->infnorm_component(iE2);
  if (debug) 
    printf("    E2 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   ENorm1,ENorm1*E2Units,EMax1,EMax1*E2Units);

  //       destroy HYPRE solver and preconditioner structures
  HYPRE_StructPCGDestroy(solver2);
  HYPRE_StructPFMGDestroy(preconditioner2);

  //       check for solver failure
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
    // if the final residual is too large, or is nan, quit
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {
      fprintf(stderr,"   Error: could not achieve prescribed tolerance!\n");

      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P2.mat\n");
      HYPRE_StructMatrixPrint("P2.mat",P,0);
      if (debug)  printf("Writing out rhs to file b2.vec\n");
      HYPRE_StructVectorPrint("b2.vec",rhsvec,0);
      if (debug)  printf("Writing out current solution to file x2.vec\n");
      HYPRE_StructVectorPrint("x2.vec",solvec,0);

      // dump module parameters to disk
      this->Dump(sol);

      fprintf(stderr,"  =====================================================================\n");
      ENZO_FAIL("MFSplit LinearSolve: HYPRE solver error (E2)");
    }
  }
  if (debug)  
    printf("  ---------------------------------------------------------------------\n");





  /////////// E3 equation ///////////
  //       set matrix values over grid
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, 
			3, sol, src, thisdt) != SUCCESS)
    ENZO_FAIL("MFSplit LinearSolve: SetupSystem failure (E3)");
  ilower[0] = SolvIndices[0][0];  iupper[0] = SolvIndices[0][1];
  ilower[1] = SolvIndices[1][0];  iupper[1] = SolvIndices[1][1];
  ilower[2] = SolvIndices[2][0];  iupper[2] = SolvIndices[2][1];
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 

  //       assemble matrix
  HYPRE_StructMatrixAssemble(P);

  //       insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);

  //       set the linear solver tolerance (rescale to relative residual and not actual)
  delta = sol_tolerance/rhsnorm;
  delta = min(delta, 1.0e-6);

  //       insert zero initial guess into HYPRE vector solvec
  for (ix=0; ix<size; ix++)  matentries[ix] = 0.0;
  HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, matentries);

  //       assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

//   if (debug)  printf("Writing out matrix to file P3.mat\n");
//   HYPRE_StructMatrixPrint("P3.mat",P,0);

//   if (debug)  printf("Writing out rhs to file b3.vec\n");
//   HYPRE_StructVectorPrint("b3.vec",rhsvec,0);

  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_StructSolver solver3;
  HYPRE_StructSolver preconditioner3;
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver3);
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner3);

  //          set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner3, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner3, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner3, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner3, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner3, sol_npost);

  //          set solver options
  HYPRE_StructPCGSetPrintLevel(solver3, sol_printl);
  HYPRE_StructPCGSetLogging(solver3, sol_log);
  HYPRE_StructPCGSetRelChange(solver3, 1);
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver3, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver3, 
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
		      preconditioner3);
  }
  else {    // ignore pfmg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver3, sol_maxit*500);
  }
  if (delta != 0.0)  HYPRE_StructPCGSetTol(solver3, Eflt64(delta));
  HYPRE_StructPCGSetup(solver3, P, rhsvec, solvec);

  //       solve the linear system
  if (debug)  printf("  Monochromatic solve for E3:\n");
  HYPRE_StructPCGSolve(solver3, P, rhsvec, solvec);

//   if (debug)  printf("Writing out solution vector to file x3.vec\n");
//   HYPRE_StructVectorPrint("x3.vec",solvec,0);

  //       extract solver & preconditioner statistics
  finalresid=1.0;
  Sits=0;
  Pits=0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver3, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver3, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner3, &Pits);
  totIters += Sits;
  if (debug)
    printf("    lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
	   finalresid*rhsnorm, sol_tolerance, rhsnorm, Sits, Pits);

  //       extract values from solution vector 
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  s_E = sol->GetData(iE3);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
    }
  }
  ENorm1 = sol->rmsnorm_component(iE3);
  EMax1  = sol->infnorm_component(iE3);
  if (debug) 
    printf("    E3 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   ENorm1,ENorm1*E3Units,EMax1,EMax1*E3Units);

  //       destroy HYPRE solver and preconditioner structures
  HYPRE_StructPCGDestroy(solver3);
  HYPRE_StructPFMGDestroy(preconditioner3);

  //       check for solver failure
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
    // if the final residual is too large, or is nan, quit
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {
      fprintf(stderr,"   Error: could not achieve prescribed tolerance!\n");

      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P3.mat\n");
      HYPRE_StructMatrixPrint("P3.mat",P,0);
      if (debug)  printf("Writing out rhs to file b3.vec\n");
      HYPRE_StructVectorPrint("b3.vec",rhsvec,0);
      if (debug)  printf("Writing out current solution to file x3.vec\n");
      HYPRE_StructVectorPrint("x3.vec",solvec,0);

      // dump module parameters to disk
      this->Dump(sol);

      fprintf(stderr,"  =====================================================================\n");
      ENZO_FAIL("MFSplit LinearSolve: HYPRE solver error (E3)");
    }
  }
  if (debug)  
    printf("  =====================================================================\n");

#else

  ENZO_FAIL("MFSplit Linear Solve Error: HYPRE must be enabled to use this module!");

#endif  // USE_HYPRE

  // return success
  return SUCCESS;
}
#endif  // TRANSFER
