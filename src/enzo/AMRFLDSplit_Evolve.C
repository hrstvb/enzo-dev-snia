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
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef FAST_SIB
			  SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
int SetFieldBoundaryConditions(int Field, HierarchyEntry *Grids[], 
			       int NumberOfGrids,
#ifdef FAST_SIB
			       SiblingGridList SiblingList[],
#endif
			       int level, TopGridData *MetaData,
			       ExternalBoundary *Exterior, 
			       LevelHierarchyEntry * Level);
int FindField(int f, int farray[], int n);



// This routine evolves the radiation field in an operator-split fashion, 
// subcycling the physics in the following manner: 
//     dt_rad <= dt_hydro
// Prior to completion, the routine also updates the maximum time step the 
// overall Grid module can take to meet a maximum subcycling ratio of 
// radiation to hydrodynamics.
//int AMRFLDSplit::Evolve(HierarchyEntry *ThisGrid, float dthydro)
//int AMRFLDSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, float dthydro)
int AMRFLDSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, 
			HierarchyEntry *Grids[], int NumberOfGrids,
			TopGridData *MetaData, ExternalBoundary *Exterior, 
#ifdef FAST_SIB
			SiblingGridList SiblingList[],
#endif
			float dthydro)
{

#ifdef AMR_SOLVE

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n AMRFLDSplit Evolve:\n");

  //  if (debug)  printf("AMRFLDSplit_Evolve: scaling radiation field\n");

  // scale radiation field on all relevant grids to solver units, output statistics
  LevelHierarchyEntry *Temp;
  float Etyp=0.0, Emax=0.0, dV = 1.0;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;
	  for (int dim=0; dim<rank; dim++)
	    dV *= (Temp->GridData->GetGridRightEdge(dim) -
		   Temp->GridData->GetGridLeftEdge( dim)) 
	        / n3[dim] / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	  
	  // scale radiation field
	  float *Enew = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	  for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k]/=ErScale;

	  for (int k=0; k<x0len*x1len*x2len; k++)  {
	    Etyp += Enew[k]*Enew[k]*dV;
	    Emax = max(Emax, fabs(Enew[k])*dV);
	  }

      }  // end loop over grids on this proc

  // communicate among all procs to get information on current radiation field
#ifdef USE_MPI
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg one = 1;
  float dtmp;
  MPI_Allreduce(&Emax, &dtmp, one, FDataType, MPI_MAX, MPI_COMM_WORLD);
  Emax = dtmp;
  MPI_Allreduce(&Etyp, &dtmp, one, FDataType, MPI_SUM, MPI_COMM_WORLD);
  Etyp = sqrt(dtmp);
#else
  Etyp = sqrt(Etyp);
#endif
  if (debug) {
    printf("   current internal quantities:\n");
    printf("      Eg rms = %13.7e, max = %13.7e\n", Etyp, Emax);
  }    


  //  if (debug)  printf("AMRFLDSplit_Evolve: attaching to amrsolve_hierarchy\n");

  // insert Enzo grids into an AMRsolve hierarchy
  AMRsolve_Hierarchy* hierarchy = new AMRsolve_Hierarchy;
//  hierarchy->enzo_attach_fld(LevelArray, level, MAX_DEPTH_OF_HIERARCHY);
  hierarchy->enzo_attach_fld(LevelArray, 0, MAX_DEPTH_OF_HIERARCHY);

  // Initialize the AMRsolve domain & hierarchy
  AMRsolve_Domain domain(3, DomainLeftEdge, DomainRightEdge);
  bool is_periodic[] = {BdryType[0][0]==0, BdryType[1][0]==0, BdryType[2][0]==0};
  hierarchy->initialize(domain, *pmpi, is_periodic);

  // initialize variables that we'll use throughout the time subcycling
  float stime2, ftime2;   // radiation, chemistry timers
  float thisdt;   // chemistry time-stepping variables
  Eflt64 Echange;
  int radstep, radstop;  // subcycle iterators


  ////////////////////////////////////
  // Problem Solve Phase

  //  if (debug)  printf("AMRFLDSplit_Evolve: getting current time\n");

  // Find a grid on this level, owned by this processor, to get current time
  HierarchyEntry* ThisGrid;
  for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
    ThisGrid = Temp->GridHierarchyEntry;
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
      break;
  }

  // Get general time-related information
  tnew = ThisGrid->GridData->ReturnTime();

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
      
      //  if (debug)  printf("AMRFLDSplit_Evolve: calling RadStep\n");

      // take a radiation step
      recompute_step = this->RadStep(LevelArray, level, hierarchy, &Echange);

      // if the radiation step was unsuccessful, update dtrad and try again
      if (recompute_step)  dtrad = max(dtrad/10, mindt);

    }

	
    // stop MPI timer for radiation solver, increment total
#ifdef USE_MPI
    ftime2 = MPI_Wtime();
#else
    ftime2 = 0.0;
#endif
    AMRSolTime += ftime2-stime2;
    
    //  if (debug)  printf("AMRFLDSplit_Evolve: computing new radiation time step\n");

    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est = this->ComputeTimeStep(Echange);
    dtrad = min(dt_est, 1.1*dtrad);

    //  if (debug)  printf("AMRFLDSplit_Evolve: updating neighbor information\n");

    //////////////////////////////////////////////////////////////////////
    // UPDATE THE FOLLOWING -- THIS CURRENTLY ONLY WORKS FOR THE TOPGRID!!
    // [also only seems to work in serial]

//     // have U0 communicate neighbor information
//     int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
//     int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
//     int ghXl = DEFAULT_GHOST_ZONES;
//     int n3[] = {1, 1, 1};
//     int dim;
//     for (dim=0; dim<rank; dim++)
//       n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
//  	      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
//     int x0len = n3[0] + 2*ghXl;
//     int x1len = n3[1] + 2*ghYl;
//     int x2len = n3[2] + 2*ghZl;
//     int NBors[3][2];
//     for (dim=0; dim<rank; dim++) 
//       for (int face=0; face<2; face++) 
// 	NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);
//     for (dim=0; dim<rank; dim++) {
//       if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
// 	NBors[dim][0] = MPI_PROC_NULL;
//       if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
// 	NBors[dim][1] = MPI_PROC_NULL;
//     }
//     EnzoVector* U0 = new EnzoVector(n3[0], n3[1], n3[2], ghXl, ghXl,
// 				    ghYl, ghYl, ghZl, ghZl, 1, NBors[0][0], 
// 				    NBors[0][1], NBors[1][0], NBors[1][1], 
// 				    NBors[2][0], NBors[2][1], 1);
//     U0->SetData(0, ThisGrid->GridData->AccessRadiationFrequency0());
//     if (U0->exchange() != SUCCESS) 
//       ENZO_FAIL("AMRFLDSplit Evolve: vector exchange error");

    //////////////////////////////////////////////////////////////////////
    // Initial attempt at just reusing Enzo's SetBoundaryConditions 
    // routine to handle ghost zone and overlap exchanges.

    if (SetBoundaryConditions(Grids, NumberOfGrids, 
#ifdef FAST_SIB
			      SiblingList, 
#endif
			      level, MetaData,
			      Exterior, NULL) == FAIL)
      ENZO_FAIL("SetBoundaryConditions() failed!\n");

    //////////////////////////////////////////////////////////////////////
    // New attempt at a rewrite of Enzo's SetBoundaryConditions 
    // routine to handle exchange of only the RadiationField

//     int FieldTypes[MAX_NUMBER_OF_BARYON_FIELDS];
//     if (ThisGrid->GridData->ReturnFieldType(FieldTypes) == FAIL) 
//       ENZO_FAIL("ReturnFieldType() failed!\n");
//     int RadField = -1;
//     for (int i=0; i<MAX_NUMBER_OF_BARYON_FIELDS; i++) 
//       if (FieldTypes[i] == RadiationFreq0) {
// 	RadField = i;
// 	break;
//       }
//     if (RadField == -1) 
//       ENZO_FAIL("Could not find RadiationFreq0 field number!\n");
// //    printf("p%"ISYM": entering SetFieldBoundaryConditions\n",MyProcessorNumber);
//     if (SetFieldBoundaryConditions(RadField, Grids, NumberOfGrids, 
// #ifdef FAST_SIB
// 				   SiblingList, 
// #endif
// 				   level, MetaData,
// 				   Exterior, NULL) == FAIL)
//       ENZO_FAIL("SetFieldBoundaryConditions() failed!\n");
// //    printf("p%"ISYM": finished SetFieldBoundaryConditions\n",MyProcessorNumber);

    

    //////////////////////////////////////////////////////////////////////



    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop
  

  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  //  if (debug)  printf("AMRFLDSplit_Evolve: detaching from amrsolve_hierarchy\n");

  // clean up amrsolve structures
  hierarchy->enzo_detach();
  delete hierarchy;
  hierarchy = NULL;

  //  if (debug)  printf("AMRFLDSplit_Evolve: filling chemistry rates\n");

  // fill the chemistry and cooling rates
  if (this->FillRates(LevelArray, level) == FAIL)
    ENZO_FAIL("AMRFLDSplit Evolve: FillRates error");

  // update the radiation time step size at this level for next time step 
  if (dtrad != huge_number) {
    for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
      for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	   Temp=Temp->NextGridThisLevel)
	Temp->GridHierarchyEntry->GridData->SetMaxRadiationDt(dtrad*maxsubcycles);
  }


  //  if (debug)  printf("AMRFLDSplit_Evolve: rescaling radiation field\n");

  // scale radiation field on all relevant grids back to enzo units
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;
	  
	  // scale radiation field
	  float *Enew = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	  for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k]*=ErScale;

      }  // end loop over grids on this proc
    


  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("AMRFLDSplit cumulative time = %g (AMRSolve = %g)\n\n",
		     RTtime, AMRSolTime);

  return SUCCESS;

#else
  ENZO_FAIL("AMRFLDSplit_Evolve ERROR: module requires AMR_SOLVE to be set!");
  return FAIL;
#endif   // AMR_SOLVE
 
}



#ifdef AMR_SOLVE

// This routine evolves the radiation subsystem within the AMRFLDSplit module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int AMRFLDSplit::RadStep(LevelHierarchyEntry *LevelArray[], int level, 
			 AMRsolve_Hierarchy *hierarchy, Eflt64 *Echange)
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
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute emissivity at this internal time step (if provided internally)
  int eta_set = 0;
#ifdef EMISSIVITY 
  if (StarMakerEmissivityField > 0)  eta_set = 1;
#endif
  if (eta_set == 0) {
    float srcNorm = this->RadiationSource(LevelArray, level, tnew);
    if (debug)   printf("   emissivity norm = %g\n",srcNorm);
  }
    
  //   enforce boundary conditions on old/new radiation fields
  if (this->EnforceBoundary(LevelArray) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   copy current radiation field into temporary field (KPhHI)
  //   on all grids owned by this processor, all levels here and down
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;

	  // access old/new radiation fields (old stored in KPhHI)
	  float *Eold = Temp->GridHierarchyEntry->GridData->AccessKPhHI();
	  float *Enew = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	  
	  // copy new radiation field into into old field
	  for (int k=0; k<x0len*x1len*x2len; k++)  Eold[k] = Enew[k];

      }  // end loop over grids on this proc
    
  // set up and solve radiation equation via amrsolve

  // Initialize the amrsolve FLD solver
  AMRsolve_Hypre_FLD amrfldsolve(*hierarchy, *amrsolve_params);
  amrfldsolve.init_hierarchy(*pmpi);
  amrfldsolve.init_stencil();
  amrfldsolve.init_graph();
//  hierarchy->print();

  //    initialize amrfldsolve system
  Eflt64 HIconst   = intSigESigHI   / intSigE;
  Eflt64 HeIconst  = intSigESigHeI  / intSigE / 4.0;
  Eflt64 HeIIconst = intSigESigHeII / intSigE / 4.0;
  Eflt64 a_ = a;
  Eflt64 a0_ = a0;
  Eflt64 adot_  = (ESpectrum == -1) ? 0.0 : adot;
  Eflt64 adot0_ = (ESpectrum == -1) ? 0.0 : adot0;
  amrfldsolve.init_elements(dt, Nchem, theta, a_, a0_, adot_, adot0_, HIconst, 
			    HeIconst, HeIIconst, NiUnits, NiUnits0, 
			    LenUnits, LenUnits0, ErUnits, ErUnits0, BdryType);

//   ENZO_FAIL("Stopping computation for debugging purposes");
  

  if (debug)  printf(" ----------------------------------------------------------------------\n");

  //    solve amrfldsolve system
  amrfldsolve.solve();
  Eflt64 finalresid = amrfldsolve.residual();
  Eint32 Sits = amrfldsolve.iterations();
  if (debug) printf("   lin resid = %.1e (tol = %.1e), its = %i\n",
		    finalresid, sol_tolerance, Sits);
  
  // check solution
  int recompute_step = 0;
  if (amrfldsolve.evaluate() != 0) {
    if (dt > mindt*TimeUnits) {
      // allow remainder of function to complete (to reset units, etc.), 
      // but have calling routine update dt and compute step again.
      recompute_step = 1;
    }
    else {
      fprintf(stderr,"AMRFLDSplit_RadStep: could not achieve prescribed tolerance!\n");
      
      // dump amrsolve matrices, module parameters to disk
      amrfldsolve.abort_dump();
      if (debug) {
	fprintf(stderr,"Dumping AMRFLDSplit module parameters to file RTdump.params\n");
	FILE *fptr = fopen("RTdump.params", "w");
	this->WriteParameters(fptr);
	fclose(fptr);
      }
      
      ENZO_FAIL("Error in AMRFLDSplit_RadStep");
    }
  }
  if (debug)  printf(" ======================================================================\n\n");

  // increment Enzo radiation field with amrsolve solution
  amrfldsolve.update_enzo();
  

  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
//  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  epsilon = 1.0e-8;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;

	  // access new radiation fields
	  float *Enew = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	  
	  // enforce floor on new field
	  for (int k=0; k<x0len*x1len*x2len; k++)  
	    Enew[k] = max(Enew[k], epsilon);

      }  // end loop over grids on this proc


  // If this was a successful solve, compute relative change in radiation field 
  // solution over time step
  if (recompute_step == 0) {
    Eflt64 pnorm = dtnorm;
    Eflt64 atol = 0.1;
    *Echange = amrfldsolve.rdiff_norm(pnorm, atol);
  }
    
  // rescale dt, told, tnew back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

}

#endif   // AMR_SOLVE


#endif   // TRANSFER
