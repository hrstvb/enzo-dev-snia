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
/  Problem initialization routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the MF solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"
#include "CosmologyParameters.h"

// character strings
EXTERN char outfilename[];


// function prototypes
int InitializeRateData(FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);

// Problem initializer prototypes



int MFProb::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

//   if (debug)  printf("Entering MFProb::Initialize routine\n");


  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) 
    ENZO_VFAIL("MFProb Initialize ERROR: p%"ISYM" could not locate his grid\n",
	       MyProcessorNumber)

  // initialize the free-streaming solver
  FSSolve = new FSProb;
  FSSolve->Initialize(TopGrid, MetaData);

  // set rank of problem
  rank = MetaData.TopGridRank;

  // get processor layout from Grid
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

  // get neighbor information from grid
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);


  // set default module parameters
  RTtime = 0.0;         // clear wall-clock time
  Nchem  = 3;           // hydrogen + helium
  Model  = 1;           // case B recombination w/ temperature
  HFrac  = 0.76;        // hydrogen fraction of total mass
  ESpectrum = 1;        // T=1e5 blackbody spectrum
  theta  = 1.0;         // backwards euler implicit time discret.
  LimImp = 0;           // lag implicit dependence of limiter in time
  LimType = 3;          // no limiter
  rScale = 1.0;         // no radiation equation scaling
  eScale = 1.0;         // no energy equation scaling
  nScale = 1.0;         // no chemistry equation scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++)   //   periodic in each direction
      BdryType[dim][face] = 0;

  // set default solver parameters
  AnalyticChem       = 0;         // theta-method
  approx_jac         = 0;         // analytical jacobian
  initial_guess      = 0;         // previous time step
  newt_maxit         = 20;        // 20 Newton iterations
  newt_norm          = 0;         // standard RMS norm
  newt_INconst       = 1.0e-4;    // inexact Newton forcing constant
  newt_tol           = 1.0e-4;    // default nonlinear tolerance
  newt_MinLinesearch = 1.0e-12;   // minimum linesearch step length
  sol_zeroguess      = 1;         // HYPRE uses a zero initial guess
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_relch          = 0;         // HYPRE relative change stopping crit.
  sol_rlxtype        = 1;         // HYPRE relaxation type
  sol_npre           = 1;         // HYPRE num pre-smoothing steps
  sol_npost          = 1;         // HYPRE num post-smoothing steps
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level

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
  float FSRadiation = 0.0;
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
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, &Nchem);
	ret += sscanf(line, "RadHydroHFraction = %"FSYM, &HFrac);
	ret += sscanf(line, "RadHydroFSRadiation = %"FSYM, &FSRadiation);
	ret += sscanf(line, "RadHydroModel = %"ISYM, &Model);
	ret += sscanf(line, "RadHydroMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "RadHydroMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "RadHydroInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "RadHydroDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "RadHydroDtRadFac = %"FSYM, &dtfac[0]);
	ret += sscanf(line, "RadHydroDtGasFac = %"FSYM, &dtfac[1]);
	ret += sscanf(line, "RadHydroDtChemFac = %"FSYM, &dtfac[2]);
	ret += sscanf(line, "RadiationScaling = %"FSYM, &rScale);
	ret += sscanf(line, "EnergyCorrectionScaling = %"FSYM, &eScale);
	ret += sscanf(line, "ChemistryScaling = %"FSYM, &nScale);
	ret += sscanf(line, "RadHydroTheta = %"FSYM, &theta);
	ret += sscanf(line, "RadHydroImplicitLimiter = %"ISYM, &LimImp);
	ret += sscanf(line, "RadHydroLimiterType = %"ISYM, &LimType);
	ret += sscanf(line, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM, 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) 
	  ret += sscanf(line, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM,
			BdryType[1], BdryType[1]+1);
	if (rank > 2) 
	  ret += sscanf(line, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM,
			BdryType[2], BdryType[2]+1);
	ret += sscanf(line, "RadHydroAprxJacobian = %"ISYM, &approx_jac);
	ret += sscanf(line, "RadHydroInitialGuess = %"ISYM, &initial_guess);
	ret += sscanf(line, "RadHydroAnalyticChem = %"ISYM, &AnalyticChem);
	ret += sscanf(line, "RadHydroSemiImplicit = %"ISYM, &semi_implicit);
	ret += sscanf(line, "RadHydroNewtLinesearch = %"ISYM, &newt_linesearch);
	ret += sscanf(line, "RadHydroNewtIters = %"ISYM, &newt_maxit);
	ret += sscanf(line, "RadHydroNewtNorm = %"ISYM, &newt_norm);
	ret += sscanf(line, "RadHydroINConst = %"FSYM, &newt_INconst);
	ret += sscanf(line, "RadHydroNewtTolerance = %"FSYM, &newt_tol);
	ret += sscanf(line, "RadHydroMinLinesearch = %"FSYM, &newt_MinLinesearch);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	
      }  // end loop over file lines

      // if doing an ionization problem (ProblemType 460),  
      // input additional parameters 
      rewind(fptr);
      if (ProblemType == 460) {
	IonizationParms[0] = 0.0;  // set defaults
	IonizationParms[1] = 0.0;
	IonizationParms[2] = 0.0;
	IonizationParms[3] = 0.0;
	IonizationParms[4] = 0.0;
        while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	  ret = 0;
	  ret += sscanf(line, "NGammaDot = %"FSYM, &IonizationParms[0]);
	  ret += sscanf(line, "EtaRadius = %"FSYM, &IonizationParms[1]);
	  ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		&IonizationParms[2], &IonizationParms[3], &IonizationParms[4]);
        }  // end loop over file lines
        if (debug) {
	  printf("MFProb_Initialize: NGammaDot = %g\n",IonizationParms[0]);
	  printf("MFProb_Initialize: EtaRadius = %g\n",IonizationParms[1]);
	  printf("MFProb_Initialize: EtaCenter = %g %g %g\n",
		 IonizationParms[2],IonizationParms[3],IonizationParms[4]);
	}
      }  // end ProblemType IF statement

    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);

  ////////////////////////////////

  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] > 2)) {
	fprintf(stderr,"MFProb_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"MFProb_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }


  // ensure that new EBdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      EBdryVals[dim][face] = NULL;


  // set up subdomain information
  //   EdgeVals gives the location of the left/right edge of the
  //      domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    EdgeVals[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    EdgeVals[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // Nchem gives the number of chemical species (must be > 0)
  if ((Nchem != 1) && (Nchem != 3)) {
    fprintf(stderr,"MFProb Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 3\n");
    Nchem = 3;  // default is hydrogen + helium
  }

  // LimImp gives the implicitness of the radiation flux limiter (see header)
  if ((LimImp < 0) || (LimImp > 4)) {
    fprintf(stderr,"MFProb Initialize: illegal LimImp = %"ISYM"\n",
	    LimImp);
    fprintf(stderr,"   re-setting LimImp to 0\n");
    LimImp = 0;  // default is time-lagged in implicit solve
  }

  // AnalyticChem redefines the gas energy and chemistry ODE functions
  if ((AnalyticChem == 1) && (approx_jac==0)) {
    fprintf(stderr,"gFLDProblem Initialize: AnalyticChem requires approximate Jacobian\n");
    fprintf(stderr,"   re-setting approx_jac to 1\n");
    approx_jac = 1;
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt gives the time step size (initialize to zero)
  dt = 0.0;

  // a, adot give cosmological expansion & rate
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;

  // *Scaling give variable scalings for implicit solver
  if (rScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal RadiationScaling = %g\n",rScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    rScale = 1.0;  // default is no scaling
  }
  if (eScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal EnergyCorrectionScaling = %g\n",eScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    eScale = 1.0;  // default is no scaling
  }
  if (nScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal ChemistryScaling = %g\n",nScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    nScale = 1.0;  // default is no scaling
  }
  if (debug)
    printf("MFProb::Initialize p%"ISYM": rScale = %g, eScale = %g, nScale = %g\n",
	   MyProcessorNumber,rScale,eScale,nScale);

  // dtfac gives the desired percent change in values per step
  if (dtfac[0] <= 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal DtRadFac = %g\n",dtfac[0]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[0] = huge_number;  // default is no limit
  }
  if (dtfac[1] <= 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal DtGasFac = %g\n",dtfac[1]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[1] = huge_number;  // default is no limit
  }
  if (dtfac[2] <= 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal DtChemFac = %g\n",dtfac[2]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[2] = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"MFProb Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 0.0 (max pointwise norm)\n");
    dtnorm = 0.0;  // default is max norm
  }

  // Theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"MFProb Initialize: illegal theta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
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
    printf("MFProb::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM"\n",
	   MyProcessorNumber, rank, Nchem);
    printf("MFProb::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
    printf("MFProb::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);
  }

  //   for non-periodic domain, unset neighbor info.
#ifndef USE_MPI
  int MPI_PROC_NULL = -3;
#endif
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  if (debug) {
    printf("MFProb::Initialize p%"ISYM": EdgeVals = (%g:%g,%g:%g,%g:%g)\n",
	   MyProcessorNumber, EdgeVals[0][0], EdgeVals[0][1], EdgeVals[1][0],
	   EdgeVals[1][1], EdgeVals[2][0], EdgeVals[2][1]);
    printf("MFProb::Initialize p%"ISYM": OnBdry = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, int(OnBdry[0][0]), int(OnBdry[0][1]),
	   int(OnBdry[1][0]),int(OnBdry[1][1]),
	   int(OnBdry[2][0]),int(OnBdry[2][1]));
    printf("MFProb::Initialize p%"ISYM": BdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber,BdryType[0][0],BdryType[0][1],
	   BdryType[1][0],BdryType[1][1],BdryType[2][0],BdryType[2][1]);
    printf("MFProb::Initialize p%"ISYM": NBors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber,NBors[0][0],NBors[0][1],
	   NBors[1][0],NBors[1][1],NBors[2][0],NBors[2][1]);
  }
  
  // set initial time step into TopGrid
  ThisGrid->GridData->SetMaxRadiationDt(initdt);
  
  // compute global dimension information
  for (dim=0; dim<rank; dim++)
    GlobDims[dim] = MetaData.TopGridDims[dim];

  // dx gives grid cell size (comoving, normalized units)
  for (dim=0; dim<rank; dim++)
    dx[dim] = (EdgeVals[dim][1]-EdgeVals[dim][0])/LocDims[dim];

  // compute global index information for this subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (EdgeVals[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      SolvIndices[dim][0] =  (long) (fCellsLeft >= 0.0) ?
	(trunc(fCellsLeft+0.5)) : (trunc(fCellsLeft-0.5));
    }

    // add on local size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + LocDims[dim] - 1;
  }

  // adjust SolvIndices, SolvOff for Dirichlet boundary zones
  for (dim=0; dim<rank; dim++) {
    SolvOff[dim] = 0;
    if (OnBdry[dim][0]  && (BdryType[dim][0] == 1)) {
      SolvIndices[dim][0] -= 1;
      SolvOff[dim] = 1;
    }
    if (OnBdry[dim][1]  && (BdryType[dim][1] == 1))
      SolvIndices[dim][1] += 1;
  }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*DEFAULT_GHOST_ZONES;

  if (debug) {
    printf("MFProb::Initialize p%"ISYM": SolvIndices = (%i:%i,%i:%i,%i:%i)\n",
	   MyProcessorNumber, SolvIndices[0][0], SolvIndices[0][1], SolvIndices[1][0], 
	   SolvIndices[1][1], SolvIndices[2][0], SolvIndices[2][1]);
    printf("MFProb::Initialize p%"ISYM": SolvOff = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, SolvOff[0], SolvOff[1], SolvOff[2]);
    printf("MFProb::Initialize p%"ISYM": LocDims = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]);
    printf("MFProb::Initialize p%"ISYM": ArrDims = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, ArrDims[0], ArrDims[1], ArrDims[2]);
  }

  // set up vector container for previous time step (empty data)
  int xghosts = DEFAULT_GHOST_ZONES, yghosts=0, zghosts=0;
  if (rank > 1)  yghosts = DEFAULT_GHOST_ZONES;
  if (rank > 2)  zghosts = DEFAULT_GHOST_ZONES;
  U0 = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], xghosts, 
		      xghosts, yghosts, yghosts, zghosts, zghosts, 
		      4+Nchem, NBors[0][0], NBors[0][1], NBors[1][0], 
		      NBors[1][1], NBors[2][0], NBors[2][1], 1);
  GhDims[0][0] = xghosts;
  GhDims[0][1] = xghosts;
  GhDims[1][0] = yghosts;
  GhDims[1][1] = yghosts;
  GhDims[2][0] = zghosts;
  GhDims[2][1] = zghosts;

  if (debug)
    printf("MFProb::Initialize p%"ISYM": GhDims = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, GhDims[0][0], GhDims[0][1], GhDims[1][0], 
	   GhDims[1][1], GhDims[2][0], GhDims[2][1]);


  // set indices for each field into EnzoVector
  iE1 = 0;     // radiation frequency 1
  iE2 = 1;     // radiation frequency 2
  iE3 = 2;     // radiation frequency 3
  iec = 3;     // gas energy correction
  iHI = 4;     // Hydrogen I
  if (Nchem == 3) {
    iHeI  = 5;   // Helium I
    iHeII = 6;   // Helium II
  }


  // set up vectors for temporary storage and Jacobian components
  sol    = U0->clone();
  extsrc = U0->clone();
  tmp1   = U0->clone();
  tmp2   = U0->clone();
  tmp3   = U0->clone();
  L_e    = U0->clone();
  L_HI   = U0->clone();
  eCorr   = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  piHI    = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  GHI     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E1_HI = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E2_HI = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E3_HI = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E1_E1 = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E2_E2 = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L_E3_E3 = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  if (Nchem == 3) {
    L_HeI  = U0->clone();
    L_HeII = U0->clone();
    piHeI     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    piHeII    = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    GHeI      = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    GHeII     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    L_E2_HeI  = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    L_E3_HeI  = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    L_E3_HeII = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  }

  
  // set up the InexactNewton solver
  INSolve = new InexactNewtonSolver(sol);


  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("MFProb Initialize: InitializeRateData failure!");
  // un-scale rates for use within RadHydro solver (handles its own units)
  {
    float MassUnits, TempUnits;
    DenUnits=LenUnits=TempUnits=TimeUnits=VelUnits=MassUnits=aUnits=1.0;
    if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	         &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) 
      ENZO_FAIL("MFProb Initialize: GetUnits failure!");
    float mp = 1.67262171e-24;   // Mass of a proton [g]
    a = 1.0; adot = 0.0;
    if (ComovingCoordinates) {
      if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL) 
	ENZO_FAIL("MFProb Initialize: CosmologyComputeExpansionFactor failure!");
      aUnits = 1.0/(1.0 + InitialRedshift);
    }
    float tbase1 = TimeUnits;
    float xbase1 = LenUnits/(a*aUnits);
    float dbase1 = DenUnits*a*a*a*aUnits*aUnits*aUnits;
    float kunit  = (aUnits*aUnits*aUnits*mp) / (dbase1*tbase1);
    float kunit_3bdy  = kunit * (aUnits*aUnits*aUnits*mp) / dbase1;
    float coolunit = (aUnits*aUnits*aUnits*aUnits*aUnits*xbase1*xbase1*mp*mp) 
                     / (tbase1*tbase1*tbase1*dbase1);
    for (i=0; i<CoolData.NumberOfTemperatureBins; i++) {
      RateData.k1[i]      *= kunit;
      RateData.k2[i]      *= kunit;
      RateData.k3[i]      *= kunit;
      RateData.k4[i]      *= kunit;
      RateData.k5[i]      *= kunit;
      RateData.k6[i]      *= kunit;
      RateData.k7[i]      *= kunit;
      RateData.k8[i]      *= kunit;
      RateData.k9[i]      *= kunit;
      RateData.k10[i]     *= kunit;
      RateData.k11[i]     *= kunit;
      RateData.k12[i]     *= kunit;
      RateData.k13[i]     *= kunit;
      RateData.k13dd[i]   *= kunit;
      RateData.k14[i]     *= kunit;
      RateData.k15[i]     *= kunit;
      RateData.k16[i]     *= kunit;
      RateData.k17[i]     *= kunit;
      RateData.k18[i]     *= kunit;
      RateData.k19[i]     *= kunit;
      RateData.k20[i]     *= kunit;
      RateData.k21[i]     *= kunit;
      RateData.k22[i]     *= kunit_3bdy;
      RateData.k50[i]     *= kunit;
      RateData.k51[i]     *= kunit;
      RateData.k52[i]     *= kunit;
      RateData.k53[i]     *= kunit;
      RateData.k54[i]     *= kunit;
      RateData.k55[i]     *= kunit;
      RateData.k56[i]     *= kunit;
      CoolData.ceHI[i]    *= coolunit;
      CoolData.ceHeI[i]   *= coolunit;
      CoolData.ceHeII[i]  *= coolunit;
      CoolData.ciHI[i]    *= coolunit;
      CoolData.ciHeI[i]   *= coolunit;
      CoolData.ciHeIS[i]  *= coolunit;
      CoolData.ciHeII[i]  *= coolunit;
      CoolData.reHII[i]   *= coolunit;
      CoolData.reHeII1[i] *= coolunit;
      CoolData.reHeII2[i] *= coolunit;
      CoolData.reHeIII[i] *= coolunit;
      CoolData.brem[i]    *= coolunit;
    }
    CoolData.comp   *= coolunit;
    CoolData.piHI   *= coolunit;
    CoolData.piHeI  *= coolunit;
    CoolData.piHeII *= coolunit;
  }


#ifdef USE_HYPRE
  // initialize HYPRE stuff
  //    initialize the diagnostic information
  totIters = 0;

  //    set up the grid
  //       create the grid object
  HYPRE_SStructGridCreate(MPI_COMM_WORLD, rank, 1, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);

  //       set grid variables
  hypre_SStructVariable_enum vartypes[3] = {HYPRE_SSTRUCT_VARIABLE_CELL,
					    HYPRE_SSTRUCT_VARIABLE_CELL,
					    HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(grid, 0, 3, vartypes);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_SStructGridSetPeriodic(grid, 0, periodicity);
  
  //       assemble the grid
  HYPRE_SStructGridAssemble(grid);

  //   set up the stencils
  if      (rank == 1)  stSize = 3;  // 3 space, 2 other species
  else if (rank == 2)  stSize = 5;  // 5 space, 2 other species
  else                 stSize = 7;  // 7 space, 2 other species
  HYPRE_SStructStencilCreate(rank, stSize+2, &stencil1);
  HYPRE_SStructStencilCreate(rank, stSize+2, &stencil2);
  HYPRE_SStructStencilCreate(rank, stSize+2, &stencil3);

  //      set stencil entries
  Eint32 offset[3];
  Eint32 stentry1=0;
  Eint32 stentry2=0;
  Eint32 stentry3=0;
  //         dependencies to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
    HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
    HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  }
  //         dependencies to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
    HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
    HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  }
  //         dependencies to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
  HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
  HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  //         dependencies to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
  HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
  HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  //         dependencies to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
  HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
  HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  //         dependencies to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
    HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
    HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  }
  //         dependencies to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 0);
    HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 1);
    HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 2);
  }
  //         dependencies to other variables
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  //            E1 dependencies on E2 and E3
  HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 1);
  HYPRE_SStructStencilSetEntry(stencil1, stentry1++, offset, 2);
  //            E2 dependencies on E1 and E3
  HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 0);
  HYPRE_SStructStencilSetEntry(stencil2, stentry2++, offset, 2);
  //            E3 dependencies on E1 and E2
  HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 0);
  HYPRE_SStructStencilSetEntry(stencil3, stentry3++, offset, 1);

  //         set up the HYPRE graph
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
  HYPRE_SStructGraphSetObjectType(graph, HYPRE_SSTRUCT);

  //         set stencils into the graph
  HYPRE_SStructGraphSetStencil(graph, 0, 0, stencil1);
  HYPRE_SStructGraphSetStencil(graph, 0, 1, stencil2);
  HYPRE_SStructGraphSetStencil(graph, 0, 2, stencil3);

  //         assemble the graph
  HYPRE_SStructGraphAssemble(graph);

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  P11tmpvec = new Eflt64[stSize*Nx*Ny*Nz];
  P12tmpvec = new Eflt64[Nx*Ny*Nz];
  P13tmpvec = new Eflt64[Nx*Ny*Nz];
  P21tmpvec = new Eflt64[Nx*Ny*Nz];
  P22tmpvec = new Eflt64[stSize*Nx*Ny*Nz];
  P23tmpvec = new Eflt64[Nx*Ny*Nz];
  P31tmpvec = new Eflt64[Nx*Ny*Nz];
  P32tmpvec = new Eflt64[Nx*Ny*Nz];
  P33tmpvec = new Eflt64[stSize*Nx*Ny*Nz];
  r1tmpvec  = new Eflt64[Nx*Ny*Nz];
  r2tmpvec  = new Eflt64[Nx*Ny*Nz];
  r3tmpvec  = new Eflt64[Nx*Ny*Nz];
  HYPREbuff = new Eflt64[Nx];
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &P);
  HYPRE_SStructMatrixSetObjectType(P, HYPRE_SSTRUCT);
  HYPRE_SStructMatrixInitialize(P);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_SStructVectorSetObjectType(rhsvec, HYPRE_SSTRUCT);
  HYPRE_SStructVectorInitialize(rhsvec);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_SStructVectorSetObjectType(solvec, HYPRE_SSTRUCT);
  HYPRE_SStructVectorInitialize(solvec);

#else

  ENZO_FAIL("MFProb Initialize Error: HYPRE must be enabled to use this module!");

#endif  // USE_HYPRE

  //   check MG solver parameters
  if (sol_maxit < 0) {
    fprintf(stderr,"Illegal RadHydroMaxMGIters = %i. Setting to 20\n",
	    sol_maxit);
    sol_maxit = 20;
  }
  if ((sol_rlxtype<0) || (sol_rlxtype>3)) {
    fprintf(stderr,"Illegal RadHydroMGRelaxType = %i. Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 0) {
    fprintf(stderr,"Illegal RadHydroMGPreRelax = %i. Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 0) {
    fprintf(stderr,"Illegal RadHydroMGPostRelax = %i. Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }

  //   check Newton solver parameters
  if (newt_maxit < 1) {
    fprintf(stderr,"Illegal RadHydroNewtIters = %"ISYM". Setting to 20\n",
	    newt_maxit);
    newt_maxit = 20;
  }
  if ((newt_norm < 0) || (newt_norm > 7)) {
    fprintf(stderr,"Illegal RadHydroNewtNorm = %"ISYM". Setting to 0\n",
	    newt_norm);
    newt_norm = 0;
  }
  if (newt_tol < 1.0e-15) {
    fprintf(stderr,"Illegal RadHydroNewtTolerance = %g. Setting to 1e-4\n",
	    newt_tol);
    newt_tol = 1.0e-4;
  }
  if ((newt_INconst <= 0.0e0) || (newt_INconst >= 1.0)) {
    fprintf(stderr,"Illegal RadHydroINConst = %g. Setting to 1.0e-4\n",
	    newt_INconst);
    newt_INconst = 1.0e-4;
  }
  if (newt_MinLinesearch < 1.0e-15) {
    fprintf(stderr,"Illegal RadHydroMinLinesearch = %g. Setting to 1e-12\n",
	    newt_MinLinesearch);
    newt_MinLinesearch = 1.0e-12;
  }

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
    
    
  // Insert new problem-specific intializers here...


  default:
    // set BC on all Dirichlet/Neumann faces to homogeneous
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x0 left radiation BC failure!");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x0 right radiation BC failure!");
    }
    if ((BdryType[1][0] != 0) && (rank > 1)) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x1 left radiation BC failure!");
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x1 right radiation BC failure!");
    }
    if ((BdryType[2][0] != 0) && (rank == 3)) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x2 left radiation BC failure!");
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFProb Initialize: x2 right radiation BC failure!");
    }
    break;
  }
  ////////////////////////////////

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      ENZO_VFAIL("Error opening parameter output file %s!\n",outfilename)
    }
    else {
      fprintf(outfptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
      fprintf(outfptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
      fprintf(outfptr, "RadHydroHFraction = %"FSYM"\n", HFrac);
      fprintf(outfptr, "RadHydroModel = %"ISYM"\n", Model);
      fprintf(outfptr, "RadHydroMaxDt = %g\n", maxdt);
      fprintf(outfptr, "RadHydroMinDt = %g\n", mindt);
      fprintf(outfptr, "RadHydroInitDt = %g\n", initdt);
      fprintf(outfptr, "RadHydroDtNorm = %"FSYM"\n", dtnorm);
      fprintf(outfptr, "RadHydroDtRadFac = %g\n", dtfac[0]);
      fprintf(outfptr, "RadHydroDtGasFac = %g\n", dtfac[1]);
      fprintf(outfptr, "RadHydroDtChemFac = %g\n", dtfac[2]);
      fprintf(outfptr, "RadiationScaling = %g\n", rScale);
      fprintf(outfptr, "EnergyCorrectionScaling = %g\n", eScale);
      fprintf(outfptr, "ChemistryScaling = %g\n", nScale);
      fprintf(outfptr, "RadHydroTheta = %g\n", theta);
      fprintf(outfptr, "RadHydroLimiterType = %"ISYM"\n", LimType);
      fprintf(outfptr, "RadHydroImplicitLimiter = %"ISYM"\n", LimImp);
      fprintf(outfptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[0][0], BdryType[0][1]);
      if (rank > 1) {
	fprintf(outfptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
		BdryType[1][0], BdryType[1][1]);
	if (rank > 2) {
	  fprintf(outfptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
		  BdryType[2][0], BdryType[2][1]);
	}
      }
      fprintf(outfptr, "RadHydroAprxJacobian = %"ISYM"\n", approx_jac);    
      fprintf(outfptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);
      fprintf(outfptr, "RadHydroAnalyticChem = %"ISYM"\n", AnalyticChem);
      fprintf(outfptr, "RadHydroSemiImplicit = %"ISYM"\n", semi_implicit);
      fprintf(outfptr, "RadHydroNewtLinesearch = %"ISYM"\n", newt_linesearch);
      fprintf(outfptr, "RadHydroNewtIters = %"ISYM"\n", newt_maxit);    
      fprintf(outfptr, "RadHydroNewtNorm = %"ISYM"\n", newt_norm);    
      fprintf(outfptr, "RadHydroINConst = %g\n", newt_INconst);    
      fprintf(outfptr, "RadHydroNewtTolerance = %g\n", newt_tol);    
      fprintf(outfptr, "RadHydroMinLinesearch = %g\n", newt_MinLinesearch);
      fprintf(outfptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "RadHydroMGPostRelax = %i\n", sol_npost);    
      
      // close parameter file
      fclose(outfptr);
    }
  }

  // set prepared flag
  prepared = true;

  return SUCCESS;
}
#endif
