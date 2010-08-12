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
#include "MFSplit.h"
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
int MFIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local,
			       float FSRadiation, float Radiation1, 
			       float Radiation2, float Radiation3,
			       float E1Units, float E2Units, 
			       float E3Units);



int MFSplit::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

//   if (debug)  printf("Entering MFSplit::Initialize routine\n");


  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) 
    ENZO_VFAIL("MFSplit Initialize ERROR: p%"ISYM" could not locate his grid\n",
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
  HFrac  = 0.76;        // hydrogen fraction of total matter
  ESpectrum = 1;        // T=1e5 blackbody spectrum
  theta  = 1.0;         // backwards euler implicit time discret.
  LimType = 3;          // no limiter
  rScale = 1.0;         // no radiation equation scaling
  eScale = 1.0;         // no energy equation scaling
  nScale = 1.0;         // no chemistry equation scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++)   //   periodic in each direction
      BdryType[dim][face] = 0;

  // set default solver parameters
  initial_guess      = 0;         // previous time step
  sol_tolerance      = 1.0e-4;    // HYPRE solver tolerance tolerance
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_rlxtype        = 1;         // HYPRE relaxation type
  sol_npre           = 1;         // HYPRE num pre-smoothing steps
  sol_npost          = 1;         // HYPRE num post-smoothing steps
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level

  // set default ionization parameters
  NGammaDot    = 0.0;             // ionization strength
  EtaRadius    = 0.0;             // single cell
  EtaCenter[0] = 0.0;             // x-location
  EtaCenter[1] = 0.0;             // y-location
  EtaCenter[2] = 0.0;             // z-location

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
	ret += sscanf(line, "RadHydroLimiterType = %"ISYM, &LimType);
	ret += sscanf(line, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM, 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) 
	  ret += sscanf(line, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM,
			BdryType[1], BdryType[1]+1);
	if (rank > 2) 
	  ret += sscanf(line, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM,
			BdryType[2], BdryType[2]+1);
	ret += sscanf(line, "RadHydroInitialGuess = %"ISYM, &initial_guess);
	ret += sscanf(line, "RadHydroNewtTolerance = %"FSYM, &sol_tolerance);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	ret += sscanf(line, "NGammaDot = %"FSYM, &NGammaDot);
	ret += sscanf(line, "EtaRadius = %"FSYM, &EtaRadius);
	ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &(EtaCenter[0]),&(EtaCenter[1]),&(EtaCenter[2]));
	
      }  // end loop over file lines
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
	fprintf(stderr,"MFSplit_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"MFSplit_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
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

  // LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // Model gives the physical set of equations to use
  if ((Model != 1) && (Model != 4)) {
    fprintf(stderr,"MFSplit Initialize: illegal Model = %"ISYM"\n",Model);
    ENZO_FAIL("   Model is unimplemented in module, exiting.");
  }

  // Nchem gives the number of chemical species (must be > 0)
  if ((Nchem != 1) && (Nchem != 3)) {
    fprintf(stderr,"MFSplit Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 3\n");
    Nchem = 3;  // default is hydrogen + helium
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt gives the time step size (initialize to zero)
  dt = 0.0;
  dtchem = huge_number;    // initially use the radiation dt

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
    printf("MFSplit::Initialize p%"ISYM": rScale = %g, eScale = %g, nScale = %g\n",
	   MyProcessorNumber,rScale,eScale,nScale);

  // dtfac gives the desired percent change in values per step
  if (dtfac[0] <= 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal DtRadFac = %g\n",dtfac[0]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[0] = huge_number;  // default is no limit
  }
  if (dtfac[1] <= 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal DtGasFac = %g\n",dtfac[1]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[1] = huge_number;  // default is no limit
  }
  if (dtfac[2] <= 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal DtChemFac = %g\n",dtfac[2]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[2] = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"MFSplit Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 0.0 (max pointwise norm)\n");
    dtnorm = 0.0;  // default is max norm
  }

  // Theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"MFSplit Initialize: illegal theta = %g\n",
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
    printf("MFSplit::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM"\n",
	   MyProcessorNumber, rank, Nchem);
    printf("MFSplit::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
    printf("MFSplit::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);
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
    printf("MFSplit::Initialize p%"ISYM": EdgeVals = (%g:%g,%g:%g,%g:%g)\n",
	   MyProcessorNumber, EdgeVals[0][0], EdgeVals[0][1], EdgeVals[1][0],
	   EdgeVals[1][1], EdgeVals[2][0], EdgeVals[2][1]);
    printf("MFSplit::Initialize p%"ISYM": OnBdry = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, int(OnBdry[0][0]), int(OnBdry[0][1]),
	   int(OnBdry[1][0]),int(OnBdry[1][1]),
	   int(OnBdry[2][0]),int(OnBdry[2][1]));
    printf("MFSplit::Initialize p%"ISYM": BdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber,BdryType[0][0],BdryType[0][1],
	   BdryType[1][0],BdryType[1][1],BdryType[2][0],BdryType[2][1]);
    printf("MFSplit::Initialize p%"ISYM": NBors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
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
    printf("MFSplit::Initialize p%"ISYM": SolvIndices = (%i:%i,%i:%i,%i:%i)\n",
	   MyProcessorNumber, SolvIndices[0][0], SolvIndices[0][1], SolvIndices[1][0], 
	   SolvIndices[1][1], SolvIndices[2][0], SolvIndices[2][1]);
    printf("MFSplit::Initialize p%"ISYM": SolvOff = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, SolvOff[0], SolvOff[1], SolvOff[2]);
    printf("MFSplit::Initialize p%"ISYM": LocDims = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]);
    printf("MFSplit::Initialize p%"ISYM": ArrDims = (%"ISYM",%"ISYM",%"ISYM")\n",
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
    printf("MFSplit::Initialize p%"ISYM": GhDims = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
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
  eCorr   = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  piHI    = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  GHI     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  if (Nchem == 3) {
    piHeI     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    piHeII    = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    GHeI      = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
    GHeII     = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  }

  
  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL)
      ENZO_FAIL("MFSplit Initialize: Error in InitializeRateData.");
  // un-scale rates for use within RadHydro solver (handles its own units)
  {
    float MassUnits, TempUnits;
    DenUnits=LenUnits=TempUnits=TimeUnits=VelUnits=MassUnits=aUnits=1.0;
    if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	         &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) {
      ENZO_FAIL("MFSplit Initialize: Error in GetUnits.");
    }
    float mp = 1.67262171e-24;   // Mass of a proton [g]
    a = 1.0; adot = 0.0;
    if (ComovingCoordinates) {
      if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error in CosmologyComputeExpansionFactor.");
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
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
  //       assemble the grid
  HYPRE_StructGridAssemble(grid);

  //   set up the stencils
  if      (rank == 1)  stSize = 3;  // 3 space, 2 other species
  else if (rank == 2)  stSize = 5;  // 5 space, 2 other species
  else                 stSize = 7;  // 7 space, 2 other species
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

  //      set stencil entries
  Eint32 offset[3];
  Eint32 stentry=0;
  //         dependencies to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependencies to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependencies to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependencies to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependencies to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependencies to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependencies to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  matentries = new Eflt64[stSize*Nx*Ny*Nz];
  rhsentries = new Eflt64[Nx*Ny*Nz];
  HYPREbuff = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

#else

  ENZO_FAIL("MFSplit Initialize Error: HYPRE must be enabled to use this module!");

#endif

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
  if (sol_tolerance < 1.0e-15) {
    fprintf(stderr,"Illegal RadHydroNewtTolerane = %g. Setting to 1e-4\n",
	    sol_tolerance);
    sol_npost = 1.0e-4;
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
    
    
  case 460:
  case 462:
    // compute compatible E1-E3 radiation values based off of FSRadiation
    float E1, E2, E3;
    if (this->RadInit(&FSRadiation, &E1, &E2, &E3, &E1Units, 
		      &E2Units, &E3Units, &chibar, &chinuint) == FAIL)
      ENZO_FAIL("MFSplit Initialize: Error in RadInit.");
    E1 = E1*1e-8;
    E2 = E2*1e-8;
    E3 = E3*1e-8;
    if (debug) {
      printf("MFSplit::Initialize: Ef=%8.2e, E1=%8.2e, E2=%8.2e, E3=%8.2e\n",
	     FSRadiation,E1,E2,E3);
      printf("MFSplit::Initialize: E*Units=%10.4e, %10.4e, %10.4e\n",
	     E1Units,E2Units,E3Units);
      printf("MFSplit::Initialize: chibar=%10.4e,  chinuint=%10.4e\n",
	     chibar,chinuint);
    }
    
    // call local problem initializer (to allocate/setup local data)
    if (MFIonizationTestInitialize(fptr, fptr, TopGrid, MetaData, 1, 
				   FSRadiation, E1, E2, E3, E1Units, 
				   E2Units, E3Units) == FAIL)
      ENZO_FAIL("MFSplit Initialize: Error in MFIonizationTestInitialize.");

    //   x0, left
    // set BC on all Dirichlet/Neumann faces to homogeneous
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x0 left radiation BCs");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x0 right radiation BCs");
    }
    if ((BdryType[1][0] != 0) && (rank > 1)) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x1 left radiation BCs");
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x1 right radiation BCs");
    }
    if ((BdryType[2][0] != 0) && (rank == 3)) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x2 left radiation BCs");
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x2 right radiation BCs");
    }
    break;
    
    
  // Insert new problem intializers here...


  default:
    // set BC on all Dirichlet/Neumann faces to homogeneous
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x0 left radiation BCs");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x0 right radiation BCs");
    }
    if ((BdryType[1][0] != 0) && (rank > 1)) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x1 left radiation BCs");
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL)
	ENZO_FAIL("MFSplit Initialize: Error setting x1 right radiation BCs");
    }
    if ((BdryType[2][0] != 0) && (rank == 3)) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x2 left radiation BCs");
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("MFSplit Initialize: Error setting x2 right radiation BCs");
    }
    break;
  }
  ////////////////////////////////

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL)
      ENZO_VFAIL("MFSplit: Error opening parameter output file %s!!\n", 
		 outfilename)
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
      fprintf(outfptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);
      fprintf(outfptr, "RadHydroNewtTolerance = %g\n", sol_tolerance);    
      fprintf(outfptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "RadHydroMGPostRelax = %i\n", sol_npost);    
      fprintf(outfptr, "NGammaDot = %g\n", NGammaDot);    
      fprintf(outfptr, "EtaRadius = %g\n", EtaRadius);    
      fprintf(outfptr, "EtaCenter = %g  %g  %g\n", 
	      EtaCenter[0], EtaCenter[1], EtaCenter[2]);    
      
      // close parameter file
      fclose(outfptr);
    }
  }

  return SUCCESS;
}
#endif
