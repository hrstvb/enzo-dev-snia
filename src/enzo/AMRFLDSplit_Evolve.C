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
/  Split Implicit Problem Class, Evolve Routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "AMRFLDSplit.h"
#include "CosmologyParameters.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);



// This routine evolves the radiation field in an operator-split fashion, 
// subcycling the physics in the following manner: 
//     dt_rad <= dt_hydro
// Prior to completion, the routine also updates the maximum time step the 
// overall Grid module can take to meet a maximum subcycling ratio of 
// radiation to hydrodynamics.
//int AMRFLDSplit::Evolve(HierarchyEntry *ThisGrid, float dthydro)
int AMRFLDSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, float dthydro)
{

//   if (debug)  printf("Entering AMRFLDSplit::Evolve routine\n");

  // Iterate over all grids on this level
  LevelHierarchyEntry *Temp;
  HierarchyEntry* ThisGrid;
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel) {
    ThisGrid = Temp->GridHierarchyEntry;

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;

#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
	    MyProcessorNumber, int(MPI_id));
    ENZO_FAIL("Error in AMRFLDSplit_Evolve");
  }
#endif

  // in case MPI is not included
#ifndef MPI_INT
  int MPI_COMM_WORLD = 0;
#endif

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n AMRFLDSplit Evolve:\n");

  // Set pointers to each variable
  float *RadiationEnergy = NULL;
  int i;
  vx = ThisGrid->GridData->AccessVelocity1();
  if (vx == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not obtain velocity1");
  vy = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not obtain velocity2");
  vz = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not obtain velocity3");
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not obtain density");
  if (DualEnergyFormalism) {
    eh = ThisGrid->GridData->AccessGasEnergy();
    if (eh == NULL) 
      ENZO_FAIL("AMRFLDSplit Evolve: could not obtain fluid energy");
  }
  else {
    eh = ThisGrid->GridData->AccessTotalEnergy();
    if (eh == NULL) 
      ENZO_FAIL("AMRFLDSplit Evolve: could not obtain fluid energy");
  }
  RadiationEnergy = ThisGrid->GridData->AccessRadiationFrequency0();
  if (RadiationEnergy == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not obtain Radiation energy");
  // "access" all chemical species (some will be NULL); this helps for 
  // problems in which we only do Hydrogen chemistry internally to this 
  // module, but Helium species are present.
  HI = HeI = HeII = NULL;
  float *HII   = NULL;
  float *HeIII = NULL;
  float *ne    = NULL;
  HI    = ThisGrid->GridData->AccessHIDensity();
  HII   = ThisGrid->GridData->AccessHIIDensity();
  HeI   = ThisGrid->GridData->AccessHeIDensity();
  HeII  = ThisGrid->GridData->AccessHeIIDensity();
  ne    = ThisGrid->GridData->AccessElectronDensity();
  //    check that we accessed the required species for Nchem
  if (Nchem > 0) 
    if (HI == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HI density!");
  if (Nchem > 1) {
    if (HeI == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HeI density!");
    if (HeII == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HeII density!");
  }

  // Get general time-related information
  tnew = ThisGrid->GridData->ReturnTime();

  // initialize external sources to 0
  extsrc->constant(0.0);

  // access/fill radiation source array (if provided externally)
  float *RadSrc = extsrc->GetData(0);
  if (RadSrc == NULL) 
    ENZO_FAIL("AMRFLDSplit Evolve: could not access Radiaton source");
  int eta_set = 0;
#ifdef EMISSIVITY
  // if using external Emissivity field source, copy into extsrc
  if (StarMakerEmissivityField > 0) {
    // access external emissivity field 
    float *EmissivitySource = ThisGrid->GridData->AccessEmissivity0();
    if (EmissivitySource == NULL) 
      ENZO_FAIL("AMRFLDSplit Evolve: could not access emissivity field");
    // copy data
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      RadSrc[i] = EmissivitySource[i];

    eta_set = 1;
    float srcNorm = extsrc->rmsnorm_component(0);
    float srcMax  = extsrc->infnorm_component(0);
    if (debug) {
      printf("   emissivity norm = %g,  max = %g\n",srcNorm,srcMax);
    }
  }
#endif

  // attach radiation array to U0 vector
  U0->SetData(0, RadiationEnergy);

  // rescale Enzo units with input scalings to non-dimensionalize within solver
  U0->scale_component(0,1.0/ErScale);

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit Evolve: vector exchange_start error");

  // output typical/maximum values
  float UTypVals = U0->rmsnorm_component(0);
  float UMaxVals = U0->infnorm_component(0);

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit Evolve: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit Evolve: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits;
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;
  if (debug) {
    printf("   current internal (physical) quantities:\n");
    printf("      Eg rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals,UTypVals*ErUnits, UMaxVals,UMaxVals*ErUnits);
  }

  // initialize variables that we'll use throughout the time subcycling
  float stime2, ftime2;   // radiation, chemistry timers
  float thisdt;   // chemistry time-stepping variables
  int radstep, radstop;  // subcycle iterators

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit Evolve: vector exchange_end error");



  ////////////////////////////////////
  // Problem Solve Phase

  // internal time-stepping loop to catch up with Hydro time
  float end_time = tnew + dthydro;
  radstop = 0;
  for (radstep=0; radstep<=maxsubcycles*10; radstep++) {
      
    // start MPI timer for radiation solver
#ifdef USE_MPI
    stime2 = MPI_Wtime();
#else
    stime2 = 0.0;
#endif

    // update time-step information
    told = tnew;

    // keep trying time steps until radiation solver succeeds. 
    // Note: if we reach the minimum time step size, RadStep will call ENZO_FAIL
    int recompute_step = 1;
    while (recompute_step) {

      // update time-step information.  Note: dtrad was set on previous 
      // iteration of solver, or by user input for first iteration
      tnew = told + dtrad;
      if ((tnew - end_time)/end_time > -1.0e-14) {   // do not exceed synchronization time
	tnew = end_time;
	radstop = 1;
      }
      dt = tnew - told;
      if (debug) 
	printf("\n subcycled rad %"ISYM": dt=%7.1e, t=%7.1e (hydro dt=%7.1e, t=%7.1e)\n",
	       radstep,dt,tnew,dthydro,end_time);
      
      // take a radiation step
      recompute_step = this->RadStep(ThisGrid, eta_set);

      // if the radiation step was unsuccessful, update dtrad and try again
      if (recompute_step)  dtrad = max(dtrad/10, mindt);

    }

	
    // stop MPI timer for radiation solver, increment total
#ifdef USE_MPI
    ftime2 = MPI_Wtime();
#else
    ftime2 = 0.0;
#endif
    HYPREtime += ftime2-stime2;
    

    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est = this->ComputeTimeStep(U0,sol);
    dtrad = min(dt_est, 1.1*dtrad);

    // update Enzo radiation field with new values
    U0->copy_component(sol, 0);
    
    // have U0 communicate neighbor information
    if (U0->exchange() != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit Evolve: vector exchange error");

    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop



  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // fill the chemistry and cooling rates
  float *phHI       = ThisGrid->GridData->AccessKPhHI();
  float *phHeI      = ThisGrid->GridData->AccessKPhHeI();
  float *phHeII     = ThisGrid->GridData->AccessKPhHeII();
  float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
  float *dissH2I    = ThisGrid->GridData->AccessKDissH2I();
  this->FillRates(sol, U0, phHI, phHeI, phHeII, photogamma, dissH2I);

  // update the radiation time step size for next time step
  if (dtrad != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(dtrad*maxsubcycles);

  // scale back to Enzo units
  U0->scale_component(0,ErScale);

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative time = %g (HYPRE = %g)\n\n",
		     RTtime, HYPREtime);

  } // for Temp = ...

  // Return
  return SUCCESS;
 
}



// This routine evolves the radiation subsystem within the AMRFLDSplit module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int AMRFLDSplit::RadStep(HierarchyEntry *ThisGrid, int eta_set)
{
   
  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits0 = RadUnits*ErScale;
  NiUnits0 = (Nchem == 0) ? 1.0 : DenUnits0/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in RadiationGetUnits.");
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits;
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew, adot to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute emissivity at this internal time step (if provided internally)
  float *RadSrc = extsrc->GetData(0);
  if (eta_set == 0) {
    if (this->RadiationSource(RadSrc, &tnew) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in RadiationSource routine");
    float srcNorm = extsrc->rmsnorm_component(0);
    float srcMax  = extsrc->infnorm_component(0);
    if (debug)  printf("   emissivity norm = %g,  max = %g\n",srcNorm,srcMax);
  }
    
  // Multigrid solver: for periodic dims, only coarsen until grid no longer divisible by 2
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
  
  
  //   enforce boundary conditions on current time step vector
  if (this->EnforceBoundary(U0) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   set initial guess as old solution
  sol->copy(U0);
    
  //   compute updated opacities
  if (this->Opacity(OpacityE, &tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in Opacity routine");
    
#ifdef USE_HYPRE
  
  HYPRE_StructSolver solver;            // HYPRE solver structure
  HYPRE_StructSolver preconditioner;    // HYPRE preconditioner structure
  
  // set up and solve radiation equation
  float *RadiationEnergy = U0->GetData(0);    // old radiation energy array
  float *Eg_new = sol->GetData(0);    // updated radiation energy array
  float rhsnorm;   // used for setting HYPRE solver tolerance
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, RadiationEnergy, 
			Eg_new, OpacityE, RadSrc) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in SetupSystem routine");
  
  // assemble matrix
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};   // matrix stencil entries
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  Eflt64 delta;    // used for setting HYPRE solver tolerance
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 
  HYPRE_StructMatrixAssemble(P);
  
  // insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);
  
  // set linear solver tolerance (rescale to relative residual and not actual)
  delta = sol_tolerance / rhsnorm;
  delta = min(delta, 1.0e-6);
  
  // insert sol initial guess into HYPRE vector x 
  int xBuff, yBuff, zBuff, Zbl, Ybl, ix, iy, iz;  // mesh indexing shortcuts
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = 0.0;
      HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
    }
  }
  
  // assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);
  
  // set up the solver [PCG] and preconditioner [PFMG]
  //    create the solver & preconditioner
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);
  
  //    set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  
  //    set solver options
  HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
  HYPRE_StructPCGSetLogging(solver, sol_log);
  HYPRE_StructPCGSetRelChange(solver, 1);
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver, 
			      (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
			      (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
			      preconditioner);
  }
  else {    // ignore smg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit*500);
  }
  if (delta != 0.0)   HYPRE_StructPCGSetTol(solver, delta);
  HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
  
  // solve the linear system
  if (debug)
    printf(" ----------------------------------------------------------------------\n");
  HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);
  
  // extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;  // HYPRE solver statistics
  Eint32 Sits=0, Pits=0;  // HYPRE solver statistics
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug) printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
		    finalresid*rhsnorm, sol_tolerance, rhsnorm, Sits, Pits);
  int recompute_step = 0;
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
    // if the final residual is too large, or is nan, set return value to reduce 
    // dt and recompute step, unless we're at the minimum step size already
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {

      if (dt > mindt*TimeUnits) {
	// allow remainder of function to complete (to reset units, etc.), 
	// but have calling routine update dt and compute step again.
	recompute_step = 1;
      }
      else {
	fprintf(stderr,"AMRFLDSplit_RadStep: could not achieve prescribed tolerance!\n");
	
	// output linear system to disk
	if (debug)  printf("Writing out matrix to file P.mat\n");
	HYPRE_StructMatrixPrint("P.mat",P,0);
	if (debug)  printf("Writing out rhs to file b.vec\n");
	HYPRE_StructVectorPrint("b.vec",rhsvec,0);
	if (debug)  printf("Writing out current solution to file x.vec\n");
	HYPRE_StructVectorPrint("x.vec",solvec,0);
	
	// dump module parameters to disk
	this->Dump(sol);
	ENZO_FAIL("Error in AMRFLDSplit_RadStep");
      }
    }
  }
  if (debug)  printf(" ======================================================================\n\n");
  
  
  // extract values from solution vector, adding them to solution
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	Eg_new[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
    }
  }
  
  // destroy HYPRE solver & preconditioner structures
  HYPRE_StructPCGDestroy(solver);
  HYPRE_StructPFMGDestroy(preconditioner);
  
#else
  ENZO_FAIL("AMRFLDSplit_RadStep ERROR: module requires USE_HYPRE to be set!");
#endif
    
  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  for (int i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    Eg_new[i] = max(Eg_new[i],epsilon);

  // rescale dt, told, tnew, adot back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;
}




#endif   // TRANSFER
