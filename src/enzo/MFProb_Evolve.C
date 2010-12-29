/***********************************************************************
/
/  Multi-Frequency, Multi-species, Implicit Problem Class
/  Evolutiona routine
/
/  written by: Daniel Reynolds
/  date:       June 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "MFProb.h"
#include "InexactNewton.h"
#include "CosmologyParameters.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);




//int MFProb::Evolve(HierarchyEntry *ThisGrid, float deltat)
int MFProb::Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat)
{

//   if (debug)  printf("Entering MFProb::Evolve routine\n");

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
  if (MyProcessorNumber != MPI_id) 
    ENZO_VFAIL("MFProb Evolve: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
	       MyProcessorNumber, int(MPI_id))
#endif

  // start MPI timer
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  // get information from Grid and inputs
  dt = deltat;
  told = ThisGrid->GridData->ReturnTime();
  tnew = told+dt;

  // Set pointers to each variable (zero out fluid energy correction 
  // since that is not internal to Enzo)
  float *Radiation1 = NULL;
  float *Radiation2 = NULL;
  float *Radiation3 = NULL;
  int i;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    eCorr[i] = 0.0;
  vx  = ThisGrid->GridData->AccessVelocity1();
  if (vx == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access velocity1!\n");
  vy  = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access velocity2!\n");
  vz  = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access velocity3!\n");
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access density!\n");
  if (DualEnergyFormalism) {
    eh  = ThisGrid->GridData->AccessGasEnergy();
    if (eh == NULL) 
      ENZO_FAIL("MFProb Evolve: could not access gas energy!\n");
  }
  else {
    eh  = ThisGrid->GridData->AccessTotalEnergy();
    if (eh == NULL) 
      ENZO_FAIL("MFProb Evolve: could not access gas energy!\n");
  }
  Efree = ThisGrid->GridData->AccessRadiationFrequency0();
  if (Efree == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access Radiation0!\n");
  Radiation1 = ThisGrid->GridData->AccessRadiationFrequency1();
  if (Radiation1 == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access Radiation1!\n");
  Radiation2 = ThisGrid->GridData->AccessRadiationFrequency2();
  if (Radiation2 == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access Radiation2!\n");
  Radiation3 = ThisGrid->GridData->AccessRadiationFrequency3();
  if (Radiation3 == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access Radiation3!\n");
  // "access" all chemical species (some will be NULL); this helps for 
  // problems in which we only do Hydrogen chemistry internally to this 
  // module, but Helium species are present.
  float *nHI    = NULL;
  float *nHII   = NULL;
  float *nHeI   = NULL;
  float *nHeII  = NULL;
  float *nHeIII = NULL;
  float *ne     = NULL;
  nHI    = ThisGrid->GridData->AccessHIDensity();
  nHII   = ThisGrid->GridData->AccessHIIDensity();
  nHeI   = ThisGrid->GridData->AccessHeIDensity();
  nHeII  = ThisGrid->GridData->AccessHeIIDensity();
  nHeIII = ThisGrid->GridData->AccessHeIIIDensity();
  ne     = ThisGrid->GridData->AccessElectronDensity();
  // we require at least Hydrogen (Nchem >= 1)
  if (nHI == NULL) 
    ENZO_FAIL("MFProb Evolve: could not access HI density!\n");
  if (Nchem == 3) {  // check for Helium species as well
    if (nHeI == NULL) 
      ENZO_FAIL("MFProb Evolve: could not access HeI density!\n");
    if (nHeII == NULL) 
      ENZO_FAIL("MFProb Evolve: could not access HeII density!\n");
  }

  // initialize external sources to 0
  extsrc->constant(0.0);

  // access/fill source arrays
  int eta_set = 0;
  float FS_NGammaDot=0.0;
#ifdef EMISSIVITY
  // If using external Emissivity field source, we will need to fill each 
  // source array with the appropriate values.  For now, Geoffrey's 
  // StarMakerEmissivityField only corresponds to the total ionizing radiation, 
  // but this will instead need to fill 4 separate arrays: 
  //    (a) one for the free-streaming portion (integrated over all frequencies)
  //    (b) one for the E1 radiation (monochromatic, at the HI ionization threshold)
  //    (c) one for the E2 radiation (monochromatic, at the HeI ionization threshold)
  //    (d) one for the E3 radiation (monochromatic, at the HeII ionization threshold)
  // Just skip this for now and issue an error, since I don't know how Geoffrey will 
  // implement his emissivity sources.
  if (StarMakerEmissivityField > 0) 
    ENZO_FAIL("MFProb Evolve: multi-frequency StarMakerEmissivity not implemented!\n");
#endif
  //   compute emissivities (if not yet set)
  if (eta_set == 0) {
    if (this->Sources(extsrc, &tnew, &FS_NGammaDot) != SUCCESS) 
      ENZO_FAIL("MFProb Evolve: Sources failure\n");
  }
  float srcNorm1 = extsrc->rmsnorm_component(iE1);
  float srcMax1  = extsrc->infnorm_component(iE1);
  float srcNorm2 = extsrc->rmsnorm_component(iE2);
  float srcMax2  = extsrc->infnorm_component(iE2);
  float srcNorm3 = extsrc->rmsnorm_component(iE3);
  float srcMax3  = extsrc->infnorm_component(iE3);
  if (debug) {
    printf("\n  MFProb Evolve emissivities (norm,max):\n");
    printf("      E1 (%7.1e,%7.1e), E2 (%7.1e,%7.1e), E3 (%7.1e,%7.1e)\n",
	   srcNorm1,srcMax1,srcNorm2,srcMax2,srcNorm3,srcMax3);
  }

  // get internal Enzo units
  //    (old time step)
  float TempUnits, MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  float RadUnits;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, told) == FAIL) 
    ENZO_FAIL("MFProb Evolve: GetUnits failure\n");
  if (RadiationGetUnits(&RadUnits, told) == FAIL) 
    ENZO_FAIL("MFProb Evolve: RadiationGetUnits failure\n");
  //    incorporate Enzo units with implicit solver unit scaling
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  rUnits0 = RadUnits*rScale/EiScale;
  eUnits  = VelUnits*VelUnits*eScale;
  nUnits0 = DenUnits0/mp*nScale;
  //    (new time step)
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, tnew) == FAIL) 
    ENZO_FAIL("MFProb Evolve: GetUnits failure\n");
  if (RadiationGetUnits(&RadUnits, tnew) == FAIL) 
    ENZO_FAIL("MFProb Evolve: RadiationGetUnits failure\n");
  // incorporate Enzo units with implicit solver unit scaling
  fsUnits = RadUnits;
  rUnits = RadUnits*rScale/EiScale;
  nUnits = DenUnits/mp*nScale;

  // set a, adot, unit scalings to correct time-level values
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) == FAIL) 
      ENZO_FAIL("MFProb Evolve: CosmologyComputeExpansionFactor failure\n");
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) == FAIL) 
      ENZO_FAIL("MFProb Evolve: CosmologyComputeExpansionFactor failure\n");
    aUnits = 1.0/(1.0 + InitialRedshift);
  }

  // adjust chemical species HI, HeI, HeII to correspond to rho (since 
  // rho and ni are advected differently in hydro solver for some reason)
  float rhochem;
  int j,k,idx;
  int flag0, flag1, flag2, flag3; flag0=flag1=flag2=flag3=0;
  // Hydrogen chemistry 
  if (HFrac != 0.0) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      // first ensure that no densities are negative
      nHI[i]  = max(0.0, nHI[i]);
      nHII[i] = max(0.0, nHII[i]);
      
      // set rhochem as the total H 'density' according to H* species
      rhochem = nHI[i] + nHII[i];
      
      // update HI as appropriate fraction of 'true' density
      nHI[i] *= rho[i]*HFrac/rhochem;
      
      // correct if precision isn't good enough, and all rho is HI
      nHI[i]  = min(nHI[i],rho[i]*HFrac);
      nHII[i] = max(0.0, rho[i]*HFrac - nHI[i]);
    }
  } else {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      nHI[i]  = 0.0;
      nHII[i] = 0.0;
    }
  }

  // Helium chemistry 
  // (or if we are running Hydrogen chemistry, but Helium is present in code)
  if ((Nchem > 1) || (HFrac < 1.0)) {
    if (HFrac != 1.0) {
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
	
	// first ensure that no densities are negative
	nHeI[i]   = max(0.0, nHeI[i]);
	nHeII[i]  = max(0.0, nHeII[i]);
	nHeIII[i] = max(0.0, nHeIII[i]);
	
	// set rhochem as the total He 'density' according to He* species
	rhochem = nHeI[i] + nHeII[i] + nHeIII[i];
	
	// update HeI, HeII as appropriate fractions of 'true' density
	nHeI[i]   *= rho[i]*(1.0-HFrac)/rhochem;
	nHeII[i]  *= rho[i]*(1.0-HFrac)/rhochem;
	nHeIII[i] *= max(0.0, rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i]);
      }
    } else {
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
	nHeI[i]   = 0.0;
	nHeII[i]  = 0.0;
	nHeIII[i] = 0.0;
      }
    }
  }

  // attach arrays to U0 vector
  U0->SetData(iE1, Radiation1);
  U0->SetData(iE2, Radiation2);
  U0->SetData(iE3, Radiation3);
  U0->SetData(iec, eCorr);
  U0->SetData(iHI, nHI);
  if (Nchem == 3) {
    U0->SetData(iHeI,  nHeI);
    U0->SetData(iHeII, nHeII);
  }

  // rescale Enzo units with input scalings to non-dimensionalize within solver
  U0->scale_component(iE1, 1.0/rScale);
  U0->scale_component(iE2, 1.0/rScale);
  U0->scale_component(iE3, 1.0/rScale);
  U0->scale_component(iec, 1.0/eScale);
  U0->scale_component(iHI, 1.0/nScale);
  if (Nchem == 3) {
    U0->scale_component(iHeI,  1.0/nScale);
    U0->scale_component(iHeII, 1.0/nScale);
  }

  // adjust monochromatic radiation densities to be bounded by 
  // free-streaming radaition (since monochromatic are advected, 
  // and free-streaming is not)
  //   enforce bounds on radiation fields
  if (this->EnforceRadiationBounds(Efree,Radiation1,Radiation2,Radiation3) == FAIL) 
    ENZO_FAIL("MFProb Evolve: EnforceRadiationBounds failure\n");

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() == FAIL) 
    ENZO_FAIL("MFProb Evolve: EnzoVector exchange_start failure\n");

  // output typical/maximum values
  float UTypVals[Nchem+4];
  float UMaxVals[Nchem+4];
  UTypVals[0] = U0->rmsnorm_component(iE1);
  UMaxVals[0] = U0->infnorm_component(iE1);
  UTypVals[1] = U0->rmsnorm_component(iE2);
  UMaxVals[1] = U0->infnorm_component(iE2);
  UTypVals[2] = U0->rmsnorm_component(iE3);
  UMaxVals[2] = U0->infnorm_component(iE3);
  UTypVals[4] = U0->rmsnorm_component(iHI);
  UMaxVals[4] = U0->infnorm_component(iHI);
  if (Nchem == 3) {
    UTypVals[5] = U0->rmsnorm_component(iHeI);
    UMaxVals[5] = U0->infnorm_component(iHeI);
    UTypVals[6] = U0->rmsnorm_component(iHeII);
    UMaxVals[6] = U0->infnorm_component(iHeII);
  }    

  //    set fluid energy correction "typical" value (since ec0=0)
  UTypVals[3] = 0.0;  UMaxVals[3] = 0.0;
  float dtmp;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
    UTypVals[3] += eh[i]*eh[i];
    UMaxVals[3] = max(UMaxVals[3],eh[i]*eh[i]);
  }
  UTypVals[3] = sqrt(UTypVals[3]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
  UMaxVals[3] = sqrt(UMaxVals[3]);
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg one = 1;
  MPI_Allreduce(&(UMaxVals[3]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  UMaxVals[3] = dtmp;
  MPI_Allreduce(&(UTypVals[3]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  UTypVals[3] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
#endif
  UTypVals[3] /= eScale;
  UMaxVals[3] /= eScale;

  if (debug) {
    printf("  MFProb Evolve: current internal (physical) quantities:\n");
    printf("      E1 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[0],UTypVals[0]*rUnits, UMaxVals[0],UMaxVals[0]*rUnits);
    printf("      E2 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[1],UTypVals[1]*rUnits, UMaxVals[1],UMaxVals[1]*rUnits);
    printf("      E3 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[2],UTypVals[2]*rUnits, UMaxVals[2],UMaxVals[2]*rUnits);
    printf("      ec rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[3],UTypVals[3]*eUnits, UMaxVals[3],UMaxVals[3]*eUnits);
    printf("      HI rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[4],UTypVals[4]*nUnits, UMaxVals[4],UMaxVals[4]*nUnits);
    if (Nchem == 3) {
      printf("     HeI rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	     UTypVals[5],UTypVals[5]*nUnits, UMaxVals[5],UMaxVals[5]*nUnits);
      printf("    HeII rms = %10.4e (%8.2e), max = %10.4e (%8.2e)",
	     UTypVals[6],UTypVals[6]*nUnits, UMaxVals[6],UMaxVals[6]*nUnits);
    }
  }

  // set prepared flag to true
  prepared = true;

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() == FAIL) 
    ENZO_FAIL("MFProb Evolve: EnzoVector exchange_end failure\n");

  // enforce boundary conditions on state U0
  if (this->EnforceBoundary(U0,0) == FAIL) 
    ENZO_FAIL("MFProb Evolve: EnforceBoundary failure\n");

  // rescale dt, told, tnew, adot to physical values
  dt *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;
  adot /= TimeUnits;
  adot0 /= TimeUnits;




  ////////////////////////////////////
  // Problem Solve Phase

  //   set the appropriate free-streaming emissivity source
  //   (unused if StarMaker emissivity is set)
  FSSolve->SetEmissivity(FS_NGammaDot);

  //   call the free-streaming solver Evolve routine
//  FSSolve->Evolve(ThisGrid, deltat);
  FSSolve->Evolve(LevelArray, level, deltat);

  //   obtain initial guess for time-evolved solution to the coupled system
  if (this->InitialGuess(sol) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("MFProb Evolve: InitialGuess failure\n");
  }

  //   set nonlinear solver parameters
  INSolve->SetDampedNewton(newt_linesearch);
  INSolve->SetMaxIters(newt_maxit);
  INSolve->SetInexactNewton(0, newt_INconst, 1.0);
  INSolve->SetNewtonTolerance(newt_tol);
  INSolve->SetNewtonNorm(newt_norm);
  INSolve->SetMinLinesearch(newt_MinLinesearch);


  // Call nonlinear solver to compute updated time step
  if (INSolve->Solve(this,sol) == FAIL) {
    // dump problem state and return failure
    this->Dump(sol);
    ENZO_FAIL("MFProb Evolve: INSolve failure\n");
  }
  // get Newton solver diagnostics
  int NewtIts = INSolve->GetNonlinearIterations();
  float FStep = INSolve->GetLinesearchStepsize();
  float FRes = INSolve->GetNonlinearResidual();
  if ((semi_implicit == 0) && (FRes > newt_tol)) {
    this->Dump(sol);
    ENZO_FAIL("MFProb Evolve: INSolve convergence failure\n");
  }


  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // update the time step size for next time step
  float RadDt = this->ComputeTimeStep(sol);
  if (RadDt != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(RadDt);

  //   enforce bounds on radiation fields
  float *sol_E1 = sol->GetData(iE1);
  float *sol_E2 = sol->GetData(iE2);
  float *sol_E3 = sol->GetData(iE3);
  if (this->EnforceRadiationBounds(Efree, sol_E1, sol_E2, sol_E3) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("MFProb Evolve: EnforceRadiationBounds failure\n");
  }

  // Rescale solution arrays to get back from solver to Enzo units
  sol->scale_component(iE1,rScale);
  sol->scale_component(iE2,rScale);
  sol->scale_component(iE3,rScale);
  sol->scale_component(iec,eScale);
  sol->scale_component(iHI,nScale);
  if (Nchem == 3) {
    sol->scale_component(iHeI, nScale);
    sol->scale_component(iHeII,nScale);
  }

  // enforce bounds on the chemical species densities
  //    Hydrogen chemistry 
  float *sol_HI = sol->GetData(iHI);
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
    sol_HI[i] = min(max(sol_HI[i],0.0),rho[i]*HFrac);
  }
  //    Helium chemistry 
  if (Nchem == 3) {
    float *sol_HeI  = sol->GetData(iHeI);
    float *sol_HeII = sol->GetData(iHeII);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      sol_HeI[i]  = min(max(0.0, sol_HeI[i]), rho[i]*(1.0-HFrac));
      sol_HeII[i] = min(max(0.0, sol_HeII[i]),rho[i]*(1.0-HFrac));
    }
  }

  // Update Enzo data with new values
  //   Radiation Energy, Chemical Species and Fluid Correction are in 
  //   sol, while Enzo pointers to these fields are in U0
  U0->copy(sol);

  //   Add fluid correction to fluid energy field (with floor)
  float *eh_tot = ThisGrid->GridData->AccessTotalEnergy();
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
    eh_tot[i] = max(eh_tot[i]+eCorr[i], tiny_number);
  if (DualEnergyFormalism) {
    float *eh_gas = ThisGrid->GridData->AccessGasEnergy();
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      eh_gas[i] = max(eh_gas[i]+eCorr[i], tiny_number);
  }

  //   Update dependent chemical species densities (ne, nHII, nHeIII) 
  //   using computed values
  if (Nchem == 1) {   // update ne, HII
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
      ne[i] = nHII[i];
    }
  }
  else if (Nchem == 3) {   // update ne, HII, HeIII
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
      nHeIII[i] = max(rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i], 0.0);
      ne[i] = nHII[i] + nHeII[i]/4.0 + nHeIII[i]/2.0;
    }
  }

  // rescale dt, told, tnew, adot to normalized values
  dt /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;
  adot *= TimeUnits;


  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative wall time = %g\n\n",RTtime);
  
  if (debug) {
    printf("  Individual timers:\n");
    for (i=0; i<30; i++)  printf("      timer %i = %g\n",i,timers[i]);
  }

  } // for Temp = ...

  // Return
  return SUCCESS;
}

#endif   // TRANSFER
