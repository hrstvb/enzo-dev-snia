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
/  Split Implicit Problem Class, Initialization routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "CosmologyParameters.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// character strings
EXTERN char outfilename[];


// function prototypes
int InitializeRateData(FLOAT Time);
int FreezeRateData(FLOAT Time, HierarchyEntry &TopGrid);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Function prototypes
int CosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid,
			      TopGridData &MetaData, int local);
int RadHydroStreamTestInitialize(FILE *fptr, FILE *Outfptr,
				 HierarchyEntry &TopGrid,
				 TopGridData &MetaData, int local);
int RadHydroRadShockInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);
int RadHydroPulseTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RadHydroConstTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);
int RHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);



int AMRFLDSplit::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

#ifdef AMR_SOLVE

  if (debug)  printf("Entering AMRFLDSplit::Initialize routine\n");

  // find root grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    ENZO_FAIL("Error in AMRFLDSplit_Initialize");
  }

#ifdef _OPENMP
  // output number of OpenMP threads that will be used in this run
  int nthreads = omp_get_max_threads();
  printf("FLD Initialize: MPI task %"ISYM" has %"ISYM" available OpenMP threads\n",
	 MyProcessorNumber,nthreads);
#endif

#ifndef MPI_INT
  // in case MPI is not included
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif

  // initialize the amrsolve Mpi object
  pmpi = new AMRsolve_Mpi(MPI_COMM_WORLD);
  AMRsolve_Grid::set_mpi(*pmpi);

  // set rank of fld problem; error message if not 3 (amrsolve requirement)
  rank = MetaData.TopGridRank;
  if (rank != 3)
    ENZO_FAIL("Error in AMRFLDSplit_Initialize: rank must be 3 (for now)");

  // initialize internal module units
  double MassUnits;
  float TempUnits;
  DenUnits=LenUnits=TempUnits=MassUnits=TimeUnits=VelUnits=aUnits=1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  a = 1.0; adot = 0.0;
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  
  // get processor layout from Grid
  int layout[3];     // number of procs in each dim (1-based)
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  int location[3];   // location of this proc in each dim (0-based)
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

//   if (debug)  printf("  Initialize: setting default parameters\n");

  // set default module parameters
  Nchem  = 1;           // hydrogen only
  int Model = 1;        // standard non-LTE, non-isothermal model
  ESpectrum = 1;        // T=10^5 blackbody spectrum
  theta  = 1.0;         // backwards euler implicit time discret.
  maxsubcycles = 1.0;   // step ratio between radiation and hydro
  dtnorm = 2.0;         // use 2-norm for time step estimation
  ErScale = 1.0;        // no radiation equation scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++)   //   periodic in each direction
      BdryType[dim][face] = 0;

  // set default solver parameters
  sol_tolerance = 1e-8;  // HYPRE solver tolerance
  sol_maxit     = 200;   // HYPRE max linear iters
  sol_type      = 1;     // HYPRE solver
  sol_printl    = 1;     // HYPRE print level
  sol_log       = 1;     // HYPRE logging level
  sol_rlxtype   = 1;     // HYPRE relaxation type
  sol_npre      = 1;     // HYPRE num pre-smoothing steps
  sol_npost     = 1;     // HYPRE num post-smoothing steps

  // set default ionization parameters
  NGammaDot    = 0.0;    // ionization strength
  EtaRadius    = 0.0;    // single cell
  EtaCenter[0] = 0.0;    // x-location
  EtaCenter[1] = 0.0;    // y-location
  EtaCenter[2] = 0.0;    // z-location

  // set default chemistry constants
  hnu0_HI   = 13.6;      // ionization energy of HI   [eV]
  hnu0_HeI  = 24.6;      // ionization energy of HeI  [eV]
  hnu0_HeII = 54.4;      // ionization energy of HeII [eV]


  ////////////////////////////////
  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ret;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "RadHydroESpectrum = %"ISYM, &ESpectrum);
	ret += sscanf(line, "RadHydroModel = %"ISYM, &Model);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, &Nchem);
	ret += sscanf(line, "RadHydroMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "RadHydroMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "RadHydroInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "RadHydroMaxSubcycles = %"FSYM, &maxsubcycles);
	ret += sscanf(line, "RadHydroDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "RadHydroDtRadFac = %"FSYM, &dtfac);
	ret += sscanf(line, "RadiationScaling = %"FSYM, &ErScale);
	ret += sscanf(line, "RadHydroTheta = %"FSYM, &theta);
	ret += sscanf(line, "RadiationBoundaryX0Faces = %i %i", 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "RadiationBoundaryX1Faces = %i %i",
			BdryType[1], BdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "RadiationBoundaryX2Faces = %i %i",
			  BdryType[2], BdryType[2]+1);
	  }
	}
	ret += sscanf(line, "RadHydroSolType = %i", &sol_type);
	ret += sscanf(line, "RadHydroSolTolerance = %"FSYM, &sol_tolerance);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	ret += sscanf(line, "NGammaDot = %"FSYM, &NGammaDot);
	ret += sscanf(line, "EtaRadius = %"FSYM, &EtaRadius);
	ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &(EtaCenter[0]), &(EtaCenter[1]), &(EtaCenter[2]));
	
      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);


  ////////////////////////////////

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;


  //// Check input parameters ////

  // First, ensure that Enzo was called with RadiativeCooling enabled 
  // (since AMRFLDSplit doesn't handle chemistry/cooling)
  if (!RadiativeCooling) {
    fprintf(stderr,"AMRFLDSplit_Initialize: RadiativeCooling must be on!  Enabling\n");
    RadiativeCooling = 1;
  }
  
  // check for appropriate BdryType values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] > 2)) {
	fprintf(stderr,"AMRFLDSplit_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"AMRFLDSplit_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }

  // ensure that new BdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      BdryVals[dim][face] = NULL;

  // Model gives the physical set of equations to use {1,4} allowed for AMRFLDSplit
  if ((Model != 1) && (Model != 4)) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal Model = %"ISYM"\n",Model);
    fprintf(stderr,"   Model is unimplemented in this module, resetting to 1.");
    Model = 1;
  }

  // Nchem gives the number of chemical species
  if ((Nchem < 1) || (Nchem > 3)) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 1\n");
    Nchem = 1;  // default is hydrogen only
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // maxsubcycles gives the maximum desired ratio between hydro time step 
  // size and radiation time step size (dt_rad <= dt_hydro)
  // ***warn if subcycling radiation***
  if (maxsubcycles < 1.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal RadHydroMaxSubcycles = %g\n",maxsubcycles);
    fprintf(stderr,"   re-setting to 1.0\n");
    maxsubcycles = 1.0;    // default is to synchronize steps
  }
  if (maxsubcycles > 1.0) {
    fprintf(stderr,"\n**************************************************************\n");
    fprintf(stderr," WARNING: radiation subcycling (RadHydroMaxSubcycles = %g > 1.0)\n",
	    maxsubcycles);
    fprintf(stderr,"          may not work properly with Enzo chemistry module!\n");
    fprintf(stderr,"**************************************************************\n\n");
  }

  // a, adot give cosmological expansion & rate
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;

  // ErScale gives variable scalings for implicit solver
  if (ErScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal RadiationScaling = %g\n",ErScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    ErScale = 1.0;  // default is no scaling
  }
  if (debug)
    printf("AMRFLDSplit::Initialize p%"ISYM": ErScale = %g\n",
	   MyProcessorNumber,ErScale);

  // dtfac gives the desired percent change in values per step
  if (dtfac <= 0.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal RadHydroDtRadFac = %g\n",dtfac);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"AMRFLDSplit Initialize: illegal theta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  //   check linear solver parameters
  if (sol_maxit < 0) {
    fprintf(stderr,"Illegal RadHydroMaxMGIters = %i. Setting to 20\n",
	    sol_maxit);
    sol_maxit = 20;
  }
  if ((sol_type < 0) || (sol_type) > 4) {
    fprintf(stderr,"Illegal RadHydroSolType = %i.  Setting to 1 (BiCGStab)\n", 
	    sol_type);
    sol_type = 1;
  }
  if ((sol_rlxtype<0) || (sol_rlxtype>3)) {
    fprintf(stderr,"Illegal RadHydroMGRelaxType = %i. Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 1) {
    fprintf(stderr,"Illegal RadHydroMGPreRelax = %i. Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 1) {
    fprintf(stderr,"Illegal RadHydroMGPostRelax = %i. Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }
  if ((sol_tolerance < 1.0e-15) || (sol_tolerance > 1.0)) {
    fprintf(stderr,"Illegal RadHydroSolTolerance = %g. Setting to 1e-4\n",
	    sol_tolerance);
    sol_tolerance = 1.0e-4;
  }


  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  if (debug){
    printf("AMRFLDSplit::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM"\n",
	   MyProcessorNumber, rank, Nchem);
    printf("AMRFLDSplit::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
  }
  if (debug) {
    printf("AMRFLDSplit::Initialize p%"ISYM": BdryType = (%i:%i,%i:%i,%i:%i)\n",
	   MyProcessorNumber, BdryType[0][0], BdryType[0][1], BdryType[1][0], 
	   BdryType[1][1], BdryType[2][0], BdryType[2][1]);
  }

  // dt* gives the time step sizes for each piece of physics
  dtrad = initdt;                        // use the input value (scaled units)

  // set a bound on the global initial dt as a factor of the radiation timestep
  dt = initdt*maxsubcycles;

  // set initial time step into TopGrid
  ThisGrid->GridData->SetMaxRadiationDt(dt);
  
  // set up vector container for previous time step (empty data)
  int xghosts = DEFAULT_GHOST_ZONES, yghosts=0, zghosts=0;
  if (rank > 1) {
    yghosts = DEFAULT_GHOST_ZONES;
    if (rank > 2) {
      zghosts = DEFAULT_GHOST_ZONES;
    }
  }

  // compute Radiation Energy spectrum integrals
  if (this->ComputeRadiationIntegrals() == FAIL) 
    ENZO_FAIL("AMRFLDSplit::Initialize Error in radiation spectrum integrals");


#ifdef USE_HYPRE

#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif
  // initialize amrsolve stuff
  //    initialize the diagnostic information
  totIters = 0;

  // set amrsolve parameters
  amrsolve_params = new AMRsolve_Parameters();
  amrsolve_params->set_defaults();
  if (sol_type == 0)  amrsolve_params->set_parameter("solver","fac");
  if (sol_type == 1)  amrsolve_params->set_parameter("solver","bicgstab");
  if (sol_type == 2)  amrsolve_params->set_parameter("solver","bicgstab-boomer");
  if (sol_type == 3)  amrsolve_params->set_parameter("solver","gmres");
  if (sol_type == 4)  amrsolve_params->set_parameter("solver","pfmg");
  char numstr[80];
  sprintf(numstr, "%e", sol_tolerance);
  amrsolve_params->set_parameter("solver_restol",numstr);
  sprintf(numstr, "%i", sol_maxit);
  amrsolve_params->set_parameter("solver_itmax",numstr);
  sprintf(numstr, "%i", sol_printl);
  amrsolve_params->set_parameter("solver_printl",numstr);
  sprintf(numstr, "%i", sol_log);
  amrsolve_params->set_parameter("solver_log",numstr);
  sprintf(numstr, "%i", sol_rlxtype);
  amrsolve_params->set_parameter("solver_rlxtype",numstr);
  sprintf(numstr, "%i", sol_npre);
  amrsolve_params->set_parameter("solver_npre",numstr);
  sprintf(numstr, "%i", sol_npost);
  amrsolve_params->set_parameter("solver_npost",numstr);
 

//   amrsolve_params->set_parameter("dump_x","true");
//   amrsolve_params->set_parameter("dump_b","true");
//   amrsolve_params->set_parameter("dump_a","true");


  if (debug) {
    printf("AMRFLDSplit::Initialize, customized amrsolve parameters:\n");
    amrsolve_params->print();
  }



#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  AMRSolTime += ftime-stime;

#else  // ifdef USE_HYPRE

  ENZO_FAIL("AMRFLDSplit_Initialize ERROR: module requires USE_HYPRE to be set!");
  
#endif

//   if (debug)  printf("  Initialize: calling local problem initializers\n");

  ////////////////////////////////
  // set up the boundary conditions on the radiation field, 
  // depending on the ProblemType
  float ZERO = 0.0;
  float ONE  = 1.0;
  float SMALL = 1.0e-6;
  fptr = NULL;

  // set boundary conditions based on problem type
  // (default to homogeneous Dirichlet)
  switch (ProblemType) {
    
  // ODE test problem, set BCs based on input.  
  // 0 implies periodic, otherwise set to homogeneous Dirichlet
  case 400:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroConstTestInitialize.");
    
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    if ( MetaData.TopGridRank >= 2 ) {
      if (BdryType[1][0] != 0) {
	if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x1 left radiation BCs.");
	if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x1 right radiation BCs.");
      }
    }
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] != 0) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Streaming test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 401:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroStreamTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroStreamTestInitialize.");
    
    //   x0, left
    if (BdryType[0][0] == 1) {
      if (this->SetupBoundary(0,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    else if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 1) {
      if (this->SetupBoundary(0,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    else if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 1) {
      if (this->SetupBoundary(1,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    else if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 1) {
      if (this->SetupBoundary(1,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    else if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] == 1) {
	if (this->SetupBoundary(2,0,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      else if (BdryType[2][0] == 2) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      
      //   x2, right
      if (BdryType[2][1] == 1) {
	if (this->SetupBoundary(2,1,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
      else if (BdryType[2][1] == 2) {
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Pulse test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 402:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroPulseTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroPulseTestInitialize.");
    
    //   x0, left
    if (BdryType[0][0] == 1) {
      if (this->SetupBoundary(0,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    else if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 1) {
      if (this->SetupBoundary(0,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    else if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 1) {
      if (this->SetupBoundary(1,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    else if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 1) {
      if (this->SetupBoundary(1,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    else if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] == 1) {
	if (this->SetupBoundary(2,0,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      else if (BdryType[2][0] == 2) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      
      //   x2, right
      if (BdryType[2][1] == 1) {
	if (this->SetupBoundary(2,1,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
      else if (BdryType[2][1] == 2) {
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Astrophysical and Lowrie & Edwards Radiating shock test problems: 
  // set Neumann value to 0.0; leave Periodic BC alone
  case 404:
  case 405:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroRadShockInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroRadShockInitialize.");
    //   x0, left
    if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if (BdryType[2][0] == 2) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[2][1] == 2) {
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 right radiation BCs.");
    }
    
    break;
    
    
  // Ionization tests 0 and 1: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 410:
  case 411:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationTestInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Ionization test 2: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 412:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationClumpInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Ionization test 13: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 413:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationSteepInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
    
  // Ionization test 14: periodic boundary conditions on all faces (store no data).
  case 414:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    break;
    


  // Ionization test 15: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 415:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Temperature jump test 16: periodic BCs on all faces
  case 416:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroConstTestInitialize.");
    break;
    

    
  // Insert new problem intializers here...


  default:
    // set BCs based on inputs, for non periodic set to 0-valued
    if (BdryType[0][0] != 0)
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    if (BdryType[0][1] != 0)
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    if (BdryType[1][0] != 0)
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    if (BdryType[1][1] != 0)
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    if (BdryType[2][0] != 0)
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 left radiation BCs.");
    if (BdryType[2][1] != 0)
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
  }
  ////////////////////////////////

  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.");
  
  // if using an isothermal "model", freeze rate data, now that ICs exist
  if (Model == 4) 
    if (FreezeRateData(MetaData.Time, TopGrid) == FAIL) 
      ENZO_FAIL("Error in FreezeRateData.");


  if (debug)  printf("  Initialize: outputting parameters to log file\n");

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      ENZO_FAIL("Error in AMRFLDSplit_Initialize");
    }
    else {
      fprintf(outfptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
      fprintf(outfptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
      fprintf(outfptr, "RadHydroMaxDt = %g\n", maxdt);
      fprintf(outfptr, "RadHydroMinDt = %g\n", mindt);
      fprintf(outfptr, "RadHydroInitDt = %g\n", initdt);
      fprintf(outfptr, "RadHydroMaxSubcycles = %g\n", maxsubcycles);
      fprintf(outfptr, "RadHydroDtNorm = %"FSYM"\n", dtnorm);
      fprintf(outfptr, "RadHydroDtRadFac = %g\n", dtfac);
      fprintf(outfptr, "RadiationScaling = %g\n", ErScale);
      fprintf(outfptr, "RadHydroTheta = %g\n", theta);
      fprintf(outfptr, "RadiationBoundaryX0Faces = %i %i\n", 
	      BdryType[0][0], BdryType[0][1]);
      if (rank > 1) {
	fprintf(outfptr, "RadiationBoundaryX1Faces = %i %i\n", 
		BdryType[1][0], BdryType[1][1]);
	if (rank > 2) {
	  fprintf(outfptr, "RadiationBoundaryX2Faces = %i %i\n", 
		  BdryType[2][0], BdryType[2][1]);
	}
      }
      fprintf(outfptr, "RadHydroSolType = %i\n", sol_type);
      fprintf(outfptr, "RadHydroSolTolerance = %g\n", sol_tolerance);    
      fprintf(outfptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "RadHydroMGPostRelax = %i\n", sol_npost);    
      fprintf(outfptr, "FSRadiationNGammaDot = %g\n", NGammaDot);
      fprintf(outfptr, "FSRadiationEtaRadius = %g\n", EtaRadius);
      fprintf(outfptr, "FSRadiationEtaCenter = %g  %g  %g\n", 
	      EtaCenter[0], EtaCenter[1], EtaCenter[2]);
      // close parameter file
      fclose(outfptr);
    }
  }

  return SUCCESS;

#else // AMR_SOLVE
  ENZO_FAIL("AMRFLDSplit_Initialize requires AMR_SOLVE in configuration");
#endif // AMR_SOLVE
}
#endif // TRANSFER
