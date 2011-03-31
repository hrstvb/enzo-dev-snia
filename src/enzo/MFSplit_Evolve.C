/***********************************************************************
/
/  Multi-Frequency, Multi-species, Split Problem Class
/  Evolution routine
/
/  written by: Daniel Reynolds
/  date:       December 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "MFSplit.h"
#include "InexactNewton.h"
#include "CosmologyParameters.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);




//int MFSplit::Evolve(HierarchyEntry *ThisGrid, float deltat)
//int MFSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat)
int MFSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, 
		    HierarchyEntry *Grids[], int NumberOfGrids,
		    TopGridData *MetaData, ExternalBoundary *Exterior, 
#ifdef FAST_SIB
		    SiblingGridList SiblingList[],
#endif
		    float deltat)
{

  // if (debug)  printf("Entering MFSplit::Evolve routine\n");

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
    ENZO_VFAIL("ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
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
    ENZO_FAIL("MFSplit Evolve: could not obtain velocity1.");
  vy  = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL)
    ENZO_FAIL("MFSplit Evolve: could not obtain velocity2.");
  vz  = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL)
    ENZO_FAIL("MFSplit Evolve: could not obtain velocity3.");
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("MFSplit Evolve: could not obtain density.");
  if (DualEnergyFormalism) {
    eh  = ThisGrid->GridData->AccessGasEnergy();
    if (eh == NULL) 
      ENZO_FAIL("MFSplit Evolve: could not obtain gas energy.");
  }
  else {
    eh  = ThisGrid->GridData->AccessTotalEnergy();
    if (eh == NULL) 
      ENZO_FAIL("MFSplit Evolve: could not obtain gas energy.");
  }
  Efree = ThisGrid->GridData->AccessRadiationFrequency0();
  if (Efree == NULL)
    ENZO_FAIL("MFSplit Evolve: could not obtain Radiation0.");
  Radiation1 = ThisGrid->GridData->AccessRadiationFrequency1();
  if (Radiation1 == NULL) 
    ENZO_FAIL("MFSplit Evolve: could not obtain Radiation1.");
  Radiation2 = ThisGrid->GridData->AccessRadiationFrequency2();
  if (Radiation2 == NULL) 
    ENZO_FAIL("MFSplit Evolve: could not obtain Radiation2.");
  Radiation3 = ThisGrid->GridData->AccessRadiationFrequency3();
  if (Radiation3 == NULL) 
    ENZO_FAIL("MFSplit Evolve: could not obtain Radiation3.");
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
    ENZO_FAIL("MFSplit Evolve: could not obtain HI density.");
  if (Nchem == 3) {  // check for Helium species as well
    if (nHeI == NULL) 
      ENZO_FAIL("MFSplit Evolve: could not obtain HeI density.");
    if (nHeII == NULL) 
      ENZO_FAIL("MFSplit Evolve: could not obtain HeII density.");
  }

  // get information from Grid and inputs
  dt = deltat;
  told = ThisGrid->GridData->ReturnTime();
  tnew = told+dt;


  // get internal Enzo units
  //    (old time step)
  float TempUnits, MassUnits, RadUnits;
  DenUnits0=LenUnits0=TempUnits=TimeUnits=VelUnits=MassUnits=1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, told) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: error in RadiationGetUnits.");
  //    incorporate Enzo units with implicit solver unit scaling
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  rUnits0 = RadUnits*rScale;
  nUnits0 = DenUnits0/mp*nScale;
  E1Units0 = E1Units0*rUnits0;  // prior to this, it only held chi1/chibar
  E2Units0 = E2Units0*rUnits0;  // prior to this, it only held chi2/chibar
  E3Units0 = E3Units0*rUnits0;  // prior to this, it only held chi3/chibar
  //    (new time step)
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, tnew) == FAIL)
    ENZO_FAIL("MFSplit Evolve: error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: error in RadiationGetUnits.");
  // incorporate Enzo units with implicit solver unit scaling
  fsUnits = RadUnits;
  rUnits = RadUnits*rScale;
  E1Units = E1Units*rUnits;  // prior to this, it only held chi1/chibar
  E2Units = E2Units*rUnits;  // prior to this, it only held chi2/chibar
  E3Units = E3Units*rUnits;  // prior to this, it only held chi3/chibar
  eUnits = VelUnits*VelUnits*eScale;
  nUnits = DenUnits/mp*nScale;

  // set a, adot, unit scalings to correct time-level values
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) == FAIL) 
      ENZO_FAIL("MFSplit Evolve: error in CosmologyComputeExpansionFactor.");
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) == FAIL) 
      ENZO_FAIL("MFSplit Evolve: error in CosmologyComputeExpansionFactor.");
  }


  // adjust chemical species HI, HeI, HeII to correspond to rho (since 
  // rho and ni are advected differently in hydro solver for some reason)
  float rhochem;
  // Hydrogen chemistry 
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
    // first ensure that no densities are negative
    nHI[i]  = max(0.0, nHI[i]);
    nHII[i] = max(0.0, nHII[i]);
    
    // set rhochem as the total H 'density' according to H* species
    rhochem = nHI[i] + nHII[i];
    
    // update HI as appropriate fraction of 'true' density
    nHI[i] *= rho[i]*HFrac/rhochem;
    
    // correct if precision isn't good enough, and all rho is HI
    nHI[i]  = min(nHI[i], rho[i]*HFrac);
    nHII[i] = max(0.0, rho[i]*HFrac - nHI[i]);
  }

  // Helium chemistry 
  // (or if we are running Hydrogen chemistry, but Helium is present in code)
  if (HFrac < 1.0) {
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

  // rescale dt, told, tnew, adot to physical values
  dt    *= TimeUnits;
  told  *= TimeUnits;
  tnew  *= TimeUnits;
  adot  /= TimeUnits;
  adot0 /= TimeUnits;

  // adjust monochromatic radiation densities to be bounded by 
  // free-streaming radiation (since monochromatic are advected, 
  // and free-streaming is not)
  //   enforce bounds on radiation fields
  if (this->EnforceRadiationBounds(Efree,Radiation1,Radiation2,Radiation3) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: error in EnforceRadiationBounds.");


  // initialize external sources to 0
  extsrc->constant(0.0);


  // access/fill radiation source arrays
  int eta_set = 0;
  float FS_NGammaDot=0.0;
#ifdef EMISSIVITY
  // If using external Emissivity field source, we will need to fill each 
  // source array with the appropriate values.  For now, Geoffrey's 
  // StarMakerEmissivityField only corresponds to the total ionizing radiation, 
  // but eventually this will instead need to fill 4 separate arrays: 
  //    (a) one for the free-streaming portion (integrated over all frequencies)
  //    (b) one for the E1 radiation (monochromatic, at the HI ionization threshold)
  //    (c) one for the E2 radiation (monochromatic, at the HeI ionization threshold)
  //    (d) one for the E3 radiation (monochromatic, at the HeII ionization threshold)
  // Just skip this for now and issue an error, since I don't know how Geoffrey will 
  // implement his emissivity sources.
  if (StarMakerEmissivityField > 0)
    ENZO_FAIL("MFSplit Evolve: multi-frequency StarMakerEmissivity not implemented!.");
#endif
  //   compute emissivities (if not yet set)
  if (eta_set == 0) {
    if (this->RadiationSource(extsrc, &tnew, &FS_NGammaDot) != SUCCESS)
      ENZO_FAIL("MFSplit Evolve: Error in RadiationSource.");
  }
  float srcNorm1 = extsrc->rmsnorm_component(iE1);
  float srcMax1  = extsrc->infnorm_component(iE1);
  float srcNorm2 = extsrc->rmsnorm_component(iE2);
  float srcMax2  = extsrc->infnorm_component(iE2);
  float srcNorm3 = extsrc->rmsnorm_component(iE3);
  float srcMax3  = extsrc->infnorm_component(iE3);
  if (debug) {
    printf("\n  MFSplit Evolve emissivities (norm,max):\n");
    printf("      E1 (%7.1e,%7.1e), E2 (%7.1e,%7.1e), E3 (%7.1e,%7.1e)\n",
	   srcNorm1,srcMax1,srcNorm2,srcMax2,srcNorm3,srcMax3);
  }

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() == FAIL) 
    ENZO_FAIL("MFSplit Evolve: EnzoVector exchange_start error.");

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
    printf("  MFSplit Evolve: current internal (physical) quantities:\n");
    printf("      E1 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[0],UTypVals[0]*E1Units, UMaxVals[0],UMaxVals[0]*E1Units);
    printf("      E2 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[1],UTypVals[1]*E2Units, UMaxVals[1],UMaxVals[1]*E2Units);
    printf("      E3 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
	   UTypVals[2],UTypVals[2]*E3Units, UMaxVals[2],UMaxVals[2]*E3Units);
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

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() == FAIL) 
    ENZO_FAIL("MFSplit Evolve: EnzoVector exchange_end error.");


  // // update the photo-ionization and photo-heating rates
  // if (this->ComputeRadiationIntegrals(U0) == FAIL) 
  //     ENZO_FAIL("MFSplit Evolve: ComputeRadiationIntegrals error.");

  // // output rms/max photo-ionization and photo-heating values
  // float Photomax1[6];
  // float Photorms1[6];
  // for (i=0; i<6; i++) {
  //   Photomax1[i] = 0.0;
  //   Photorms1[i] = 0.0;
  // }
  // for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
  //   Photorms1[0] += piHI[i]*piHI[i];
  //   Photorms1[3] += GHI[i]*GHI[i];

  //   Photomax1[0] = max(Photomax1[0],piHI[i]*piHI[i]);
  //   Photomax1[3] = max(Photomax1[3],GHI[i]*GHI[i]);
  // }
  // if (Nchem == 3) {
  //   for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
  //     Photorms1[1] += piHeI[i]*piHeI[i];
  //     Photorms1[2] += piHeII[i]*piHeII[i];
  //     Photorms1[4] += GHeI[i]*GHeI[i];
  //     Photorms1[5] += GHeII[i]*GHeII[i];
      
  //     Photomax1[1] = max(Photomax1[1],piHeI[i]*piHeI[i]);
  //     Photomax1[2] = max(Photomax1[2],piHeII[i]*piHeII[i]);
  //     Photomax1[4] = max(Photomax1[4],GHeI[i]*GHeI[i]);
  //     Photomax1[5] = max(Photomax1[5],GHeII[i]*GHeII[i]);
  //   }
  // }
  // for (i=0; i<6; i++) {
  //   Photorms1[i] = sqrt(Photorms1[i]);
  //   Photomax1[i] = sqrt(Photomax1[i]);
  // }
  // if (debug) {
  //   printf("  MFSplit Evolve: current photo-ionization and heating rates:\n");
  //   printf("      piHI   2norm = %10.4e, max = %10.4e\n", Photorms1[0], Photomax1[0]);
  //   if (Nchem == 3) {
  //     printf("      piHeI  2norm = %10.4e, max = %10.4e\n", Photorms1[1], Photomax1[1]);
  //     printf("      piHeII 2norm = %10.4e, max = %10.4e\n", Photorms1[2], Photomax1[2]);
  //   }
  //   printf("      GHI    2norm = %10.4e, max = %10.4e\n", Photorms1[3], Photomax1[3]);
  //   if (Nchem == 3) {
  //     printf("      GHeI   2norm = %10.4e, max = %10.4e\n", Photorms1[4], Photomax1[4]);
  //     printf("      GHeII  2norm = %10.4e, max = %10.4e\n", Photorms1[5], Photomax1[5]);
  //   }
  // }



  ////////////////////////////////////
  // Problem Solve Phase

  // enforce boundary conditions on state U0
  if (this->EnforceBoundary(U0) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: EnforceBoundary error.");

  //   set the appropriate free-streaming emissivity source
  //   (unused if StarMaker emissivity is set)
  FSSolve->SetEmissivity(FS_NGammaDot);

  //   call the free-streaming solver Evolve routine
//  FSSolve->Evolve(ThisGrid, deltat);
//  FSSolve->Evolve(LevelArray, level, deltat);
  FSSolve->Evolve(LevelArray, level, Grids, NumberOfGrids, MetaData, Exterior, 
#ifdef FAST_SIB
		  SiblingList,
#endif
		  deltat);

  //   obtain initial guess for time-evolved solution to the coupled system
  if (this->InitialGuess(sol) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("MFSplit Evolve: InitialGuess failure.");
  }

  //   use old time step as initial guess for radiation equations
  sol->copy_component(U0,iE1);
  sol->copy_component(U0,iE2);
  sol->copy_component(U0,iE3);

  //   enforce boundary conditions on new time step initial guess vector
  if (this->EnforceBoundary(sol) == FAIL)
    ENZO_FAIL("MFSplit Evolve: EnforceBoundary failure.");


  //   set up and solve linear systems for radiation equations, updating the
  //   radiation solutions in the process
  if (this->LinearSolve(sol, extsrc, dt) == FAIL) {
    // dump problem state and return failure
    this->Dump(sol);
    ENZO_FAIL("MFSplit Evolve: LinSolve failure.");
  }

  // mask the free-streaming radiation field to cut-off outside of 
  // "ionized" region
  float *sol_E1 = sol->GetData(iE1);
  float *sol_E2 = sol->GetData(iE2);
  float *sol_E3 = sol->GetData(iE3);
  if (this->FSRadiationMask(Efree, sol_E1, sol_E2, sol_E3) == FAIL) 
    ENZO_FAIL("MFSplit Evolve: error in FSRadiationMask.");

  //   enforce bounds on radiation fields
  if (this->EnforceRadiationBounds(Efree, sol_E1, sol_E2, sol_E3) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("MFSplit Evolve: EnforceRadiationBounds failure.");
  }

//   // output typical/maximum values
//   UTypVals[0] = sol->rmsnorm_component(iE1);
//   UMaxVals[0] = sol->infnorm_component(iE1);
//   UTypVals[1] = sol->rmsnorm_component(iE2);
//   UMaxVals[1] = sol->infnorm_component(iE2);
//   UTypVals[2] = sol->rmsnorm_component(iE3);
//   UMaxVals[2] = sol->infnorm_component(iE3);
//   if (debug) {
//     printf("  after radiation solve, radiation quantities:\n");
//     printf("      E1 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[0],UTypVals[0]*E1Units, UMaxVals[0],UMaxVals[0]*E1Units);
//     printf("      E2 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[1],UTypVals[1]*E2Units, UMaxVals[1],UMaxVals[1]*E2Units);
//     printf("      E3 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[2],UTypVals[2]*E3Units, UMaxVals[2],UMaxVals[2]*E3Units);
//   }


  // update the photo-ionization and photo-heating rates
  if (this->ComputeRadiationIntegrals(sol) == FAIL)
    ENZO_FAIL("MFSplit Evolve: ComputeRadiationIntegrals failure.");

  // // output rms/max photo-ionization and photo-heating values
  // float PhotoMax[6];
  // float PhotoRMS[6];
  // for (i=0; i<6; i++) {
  //   PhotoMax[i] = 0.0;
  //   PhotoRMS[i] = 0.0;
  // }
  // for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
  //   PhotoRMS[0] += piHI[i]*piHI[i];
  //   PhotoRMS[3] += GHI[i]*GHI[i];

  //   PhotoMax[0] = max(PhotoMax[0],piHI[i]*piHI[i]);
  //   PhotoMax[3] = max(PhotoMax[3],GHI[i]*GHI[i]);
  // }
  // if (Nchem == 3) {
  //   for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
  //     PhotoRMS[1] += piHeI[i]*piHeI[i];
  //     PhotoRMS[2] += piHeII[i]*piHeII[i];
  //     PhotoRMS[4] += GHeI[i]*GHeI[i];
  //     PhotoRMS[5] += GHeII[i]*GHeII[i];
      
  //     PhotoMax[1] = max(PhotoMax[1],piHeI[i]*piHeI[i]);
  //     PhotoMax[2] = max(PhotoMax[2],piHeII[i]*piHeII[i]);
  //     PhotoMax[4] = max(PhotoMax[4],GHeI[i]*GHeI[i]);
  //     PhotoMax[5] = max(PhotoMax[5],GHeII[i]*GHeII[i]);
  //   }
  // }
  // for (i=0; i<6; i++) {
  //   PhotoRMS[i] = sqrt(PhotoRMS[i]);
  //   PhotoMax[i] = sqrt(PhotoMax[i]);
  // }
  // if (debug) {
  //   printf("  MFSplit Evolve: current photo-ionization and heating rates:\n");
  //   printf("      piHI   2norm = %10.4e, max = %10.4e\n", PhotoRMS[0], PhotoMax[0]);
  //   if (Nchem == 3) {
  //     printf("      piHeI  2norm = %10.4e, max = %10.4e\n", PhotoRMS[1], PhotoMax[1]);
  //     printf("      piHeII 2norm = %10.4e, max = %10.4e\n", PhotoRMS[2], PhotoMax[2]);
  //   }
  //   printf("      GHI    2norm = %10.4e, max = %10.4e\n", PhotoRMS[3], PhotoMax[3]);
  //   if (Nchem == 3) {
  //     printf("      GHeI   2norm = %10.4e, max = %10.4e\n", PhotoRMS[4], PhotoMax[4]);
  //     printf("      GHeII  2norm = %10.4e, max = %10.4e\n", PhotoRMS[5], PhotoMax[5]);
  //   }
  // }

  // subcycle the chemistry and gas energy equations
  if (debug)  printf("  Subcycling chemistry:  dt_rad = %8.2e,  t_rad = %8.2e\n", 
		     dt/TimeUnits, tnew/TimeUnits);
  float thisdt, dtchem2;
  float tchem = told;
  float *ecsrc   = extsrc->GetData(iec);
  float *HIsrc   = extsrc->GetData(iHI);
  float *HeIsrc  = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  float *sol_ec   = sol->GetData(iec);
  float *sol_HI   = sol->GetData(iHI);
  float *sol_HeI  = sol->GetData(iHeI);
  float *sol_HeII = sol->GetData(iHeII);
  float *tmp, factor, epsilon2, *eh_tot, *eh_gas;
  for (int chemstep=0; chemstep<=100; chemstep++) {
    // update tchem
    thisdt = min(dtchem, dt);            // do not exceed radiation dt
    thisdt = max(thisdt, dt/100);        // take at most 100 steps
    tchem += thisdt;                     // update chemistry time
    //    check that we don't exceed radiation time
    if (tchem >= tnew) {
      thisdt = tnew - (tchem - thisdt);  // get max time step
      tchem  = tnew;                     // set updated time
    }
    if (debug)  printf("    chem step %"ISYM":  dt_chem = %8.2e, t_chem = %8.2e\n",
	     chemstep, thisdt/TimeUnits, tchem/TimeUnits);

    //    fill in the gas energy and chemistry source terms
    if (this->GasEnergySource(extsrc, &tchem) == FAIL)
      ENZO_FAIL("MFSplit Evolve: GasEnergySource failure.");
    if (this->ChemistrySource(extsrc, &tchem) == FAIL) 
      ENZO_FAIL("MFSplit Evolve: ChemistrySource failure.");
      

    //    solve local chemistry/gas energy systems
    if (this->AnalyticChemistry(sol, extsrc, thisdt) != SUCCESS) 
      ENZO_FAIL("MFSplit Evolve: AnalyticChemistry failure.");

    //    update chemistry time step size based on changes to chem+energy
    dtchem2 = this->ComputeTimeStep(sol,1)*TimeUnits;
    dtchem  = min(dtchem2, 2.0*dtchem);

    //    enforce a solution floor on number densities
    epsilon2 = 0.0;       // put a hard floor of 0 on these fields
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      sol_HI[i] = min(max(sol_HI[i],epsilon2),rho[i]*HFrac);
    if (Nchem == 3) {
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
	sol_HeI[i]  = max(sol_HeI[i],epsilon2);
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
	sol_HeII[i] = max(sol_HeII[i],epsilon2);
    }

    //    add fluid correction to fluid energy field (with floor)
    eh_tot = ThisGrid->GridData->AccessTotalEnergy();
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      eh_tot[i] = max(eh_tot[i]+sol_ec[i]*eScale, tiny_number);
    if (DualEnergyFormalism) {
      eh_gas = ThisGrid->GridData->AccessGasEnergy();
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
	eh_gas[i] = max(eh_gas[i]+sol_ec[i]*eScale, tiny_number);
    }
    
    // Update Enzo chemistry arrays with new values
    U0->copy_component(sol, iHI);
    if (Nchem == 3) {
      U0->copy_component(sol, iHeI);
      U0->copy_component(sol, iHeII);
    }

    // break out of time-stepping loop if we've reached the end
    if (tchem >= tnew)  break;

  }
  

  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // update the time step size for next time step
  float RadDt = this->ComputeTimeStep(sol,0);
  if (debug)  printf("  MFSplit time step = %g\n",RadDt);
  if (RadDt != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(RadDt);

  // Update Enzo radiation data with new values
  U0->copy_component(sol, iE1);
  U0->copy_component(sol, iE2);
  U0->copy_component(sol, iE3);



//   // output typical/maximum values
//   UTypVals[0] = sol->rmsnorm_component(iE1);
//   UMaxVals[0] = sol->infnorm_component(iE1);
//   UTypVals[1] = sol->rmsnorm_component(iE2);
//   UMaxVals[1] = sol->infnorm_component(iE2);
//   UTypVals[2] = sol->rmsnorm_component(iE3);
//   UMaxVals[2] = sol->infnorm_component(iE3);
//   UTypVals[4] = sol->rmsnorm_component(iHI);
//   UMaxVals[4] = sol->infnorm_component(iHI);
//   if (Nchem == 3) {
//     UTypVals[5] = sol->rmsnorm_component(iHeI);
//     UMaxVals[5] = sol->infnorm_component(iHeI);
//     UTypVals[6] = sol->rmsnorm_component(iHeII);
//     UMaxVals[6] = sol->infnorm_component(iHeII);
//   }    

//   //    set fluid energy correction "typical" value (since ec0=0)
//   UTypVals[3] = 0.0;  UMaxVals[3] = 0.0;
//   for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
//     UTypVals[3] += eh[i]*eh[i];
//     UMaxVals[3] = max(UMaxVals[3],eh[i]*eh[i]);
//   }
//   UTypVals[3] = sqrt(UTypVals[3]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
//   UMaxVals[3] = sqrt(UMaxVals[3]);
// #ifdef USE_MPI
//   MPI_Allreduce(&(UMaxVals[3]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
//   UMaxVals[3] = dtmp;
//   MPI_Allreduce(&(UTypVals[3]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
//   UTypVals[3] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
// #endif
//   UTypVals[3] /= eScale;
//   UMaxVals[3] /= eScale;

//   if (debug) {
//     printf("  resulting internal (physical) quantities:\n");
//     printf("      E1 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[0],UTypVals[0]*E1Units, UMaxVals[0],UMaxVals[0]*E1Units);
//     printf("      E2 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[1],UTypVals[1]*E2Units, UMaxVals[1],UMaxVals[1]*E2Units);
//     printf("      E3 rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[2],UTypVals[2]*E3Units, UMaxVals[2],UMaxVals[2]*E3Units);
//     printf("      ec rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[3],UTypVals[3]*eUnits, UMaxVals[3],UMaxVals[3]*eUnits);
//     printf("      HI rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	   UTypVals[4],UTypVals[4]*nUnits, UMaxVals[4],UMaxVals[4]*nUnits);
//     if (Nchem == 3) {
//       printf("     HeI rms = %10.4e (%8.2e), max = %10.4e (%8.2e)\n",
// 	     UTypVals[5],UTypVals[5]*nUnits, UMaxVals[5],UMaxVals[5]*nUnits);
//       printf("    HeII rms = %10.4e (%8.2e), max = %10.4e (%8.2e)",
// 	     UTypVals[6],UTypVals[6]*nUnits, UMaxVals[6],UMaxVals[6]*nUnits);
//     }
//   }



  // Rescale solution arrays to get back from solver to Enzo units
  U0->scale_component(iE1,rScale);
  U0->scale_component(iE2,rScale);
  U0->scale_component(iE3,rScale);
  U0->scale_component(iec,eScale);
  U0->scale_component(iHI,nScale);
  if (Nchem == 3) {
    U0->scale_component(iHeI, nScale);
    U0->scale_component(iHeII,nScale);
  }

  //   Update dependent chemical species densities (ne, nHII, nHeIII) 
  //   using computed values
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
    nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
  if (Nchem == 3) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      nHeIII[i] = max(rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i], 0.0);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      ne[i] = nHII[i] + nHeII[i]/4.0 + nHeIII[i]/2.0;
  }
  else 
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      ne[i] = nHII[i];

  // rescale dt, told, tnew, adot to normalized values
  dt    /= TimeUnits;
  told  /= TimeUnits;
  tnew  /= TimeUnits;
  adot  *= TimeUnits;
  adot0 *= TimeUnits;

  // rescale E*Units to remove redshift-dependent component so that 
  // it can be reused in the next time step
  E1Units = E1Units/rUnits;
  E2Units = E2Units/rUnits;
  E3Units = E3Units/rUnits;
  E1Units0 = E1Units0/rUnits0;
  E2Units0 = E2Units0/rUnits0;
  E3Units0 = E3Units0/rUnits0;

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("  MFSplit cumulative wall time = %g\n\n",RTtime);
  
  } // for Temp = ...

  // Return
  return SUCCESS;
 
}

#endif   // TRANSFER
