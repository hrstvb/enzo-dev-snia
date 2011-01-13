/***********************************************************************
/
/  INITIALIZE MULTI-FREQUENCY RAD-HYDRO IONIZATION TEST
/
/  written by: Daniel Reynolds
/  date:       September 2009
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


/* default constants */
#define MIN_TEMP 1.0             // minimum temperature [K]
#define MAX_INITIAL_PATCHES 100  // max number of parameter file defined subgrids


// function prototypes
int InitializeRateData(FLOAT Time);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);



int MFIonizationTestInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local,
			       float FSRadiation, float Radiation1, 
			       float Radiation2, float Radiation3, 
			       float E1Units, float E2Units, float E3Units)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering MFIonizationTestInitialize routine\n");

  char *kphHIName    = "HI_kph";
  char *kphHeIName   = "HeI_kph";
  char *kphHeIIName  = "HeII_kph";
  char *gammaName    = "PhotoGamma";
  char *kdissH2IName = "H2I_kdiss";
  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *FSRadName = "FS_Radiation_Energy";
  char *Rad1Name  = "Radiation_Energy1";
  char *Rad2Name  = "Radiation_Energy2";
  char *Rad3Name  = "Radiation_Energy3";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *DeName    = "Electron_Density";
  char *EtaName   = "Emissivity";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, gridnum, ret, level, patch;

  // Setup and parameters:
  //   ambient density (should be very small) - free parameter
  //   ambient gas velocity - free parameter
  //   ambient gas temperature
  //   Hydrogen mass fraction 
  //   initial fraction HII
  //   initial fraction HeII
  //   initial fraction HeIII
  //   Number of chemical species
  //   mesh spacing
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroDensity              = 10.0;
  float RadHydroTemperature          = 1.0;
  float RadHydroIEnergy              = -1.0;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  float RadHydroInitialFractionHeII  = 0.0;
  float RadHydroInitialFractionHeIII = 0.0;
  int   RadHydroChemistry            = 3;
  int   RadHydroModel                = 1;
  int   AMRNumberOfInitialPatches = 0;
  int   AMRPatchLevel[MAX_INITIAL_PATCHES];
  FLOAT AMRPatchLeftEdge[MAX_INITIAL_PATCHES][MAX_DIMENSION];
  FLOAT AMRPatchRightEdge[MAX_INITIAL_PATCHES][MAX_DIMENSION];
  for (i=0; i<MAX_INITIAL_PATCHES; i++) {
    AMRPatchLevel[i] = 0;
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      AMRPatchLeftEdge[i][dim] = 0.0;
      AMRPatchRightEdge[i][dim] = 0.0;
    }
  }

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "RadHydroVelocity = %"FSYM" %"FSYM" %"FSYM,
		      &RadHydroX0Velocity, &RadHydroX1Velocity, 
		      &RadHydroX2Velocity);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, 
		      &RadHydroChemistry);
	ret += sscanf(line, "RadHydroModel = %"ISYM, 
		      &RadHydroModel);
	ret += sscanf(line, "RadHydroDensity = %"FSYM, 
		      &RadHydroDensity);
	ret += sscanf(line, "RadHydroTemperature = %"FSYM, 
		      &RadHydroTemperature);
	ret += sscanf(line, "RadHydroIEnergy = %"FSYM, 
		      &RadHydroIEnergy);
	if (RadHydroChemistry > 0) {
	  ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
			&RadHydroInitialFractionHII);
	  ret += sscanf(line, "RadHydroHFraction = %"FSYM, 
			&RadHydroHydrogenMassFraction);
	  ret += sscanf(line, "RadHydroInitialFractionHeII = %"FSYM, 
			&RadHydroInitialFractionHeII);
	  ret += sscanf(line, "RadHydroInitialFractionHeIII = %"FSYM, 
			&RadHydroInitialFractionHeIII);
	}

	// AMR hierarchy information
	ret += sscanf(line, "AMRNumberOfInitialPatches = %"ISYM, &AMRNumberOfInitialPatches);
	if (sscanf(line, "AMRPatchLevel[%"ISYM"]", &gridnum) > 0)
	  ret += sscanf(line, "AMRPatchLevel[%"ISYM"] = %"ISYM,
			&gridnum, &AMRPatchLevel[gridnum]);
	if (sscanf(line, "AMRPatchLeftEdge[%"ISYM"]", &gridnum) > 0)
	  ret += sscanf(line, "AMRPatchLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
			&gridnum, &AMRPatchLeftEdge[gridnum][0],
			&AMRPatchLeftEdge[gridnum][1], &AMRPatchLeftEdge[gridnum][2]);
	if (sscanf(line, "AMRPatchRightEdge[%"ISYM"]", &gridnum) > 0)
	  ret += sscanf(line, "AMRPatchRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
			&gridnum, &AMRPatchRightEdge[gridnum][0],
			&AMRPatchRightEdge[gridnum][1], &AMRPatchRightEdge[gridnum][2]);
 
      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  // check AMR inputs
  if (AMRNumberOfInitialPatches > MAX_INITIAL_PATCHES) {
    ENZO_FAIL("Too many InitialPatches! increase MAX_INITIAL_PATCHES\n");
  }

  // output requested initial AMR hierarchy structure
  if (debug && local && AMRNumberOfInitialPatches) {
    printf("  Initializing with prescribed AMR hierarchy:\n");
    printf("    AMRNumberOfInitialPatches = %"ISYM"\n",
	   AMRNumberOfInitialPatches);
    for (patch=0; patch<AMRNumberOfInitialPatches; patch++) {
      printf("    AMRPatchLevel[%"ISYM"] = %"ISYM"\n", patch, 
	     AMRPatchLevel[patch]);
      printf("    AMRPatchLeftEdge[%"ISYM"] = %g %g %g\n", patch, 
	     AMRPatchLeftEdge[patch][0], AMRPatchLeftEdge[patch][1], 
	     AMRPatchLeftEdge[patch][2]);
      printf("    AMRPatchRightEdge[%"ISYM"] = %g %g %g\n", patch, 
	     AMRPatchRightEdge[patch][0], AMRPatchRightEdge[patch][1], 
	     AMRPatchRightEdge[patch][2]);
    }
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("MFIonizationTestInitialize: InitializeRateData failure!\n");

  // if temperature specified and not internal energy, perform conversion here
  if (RadHydroIEnergy == -1.0) {
    if (RadHydroTemperature == -1.0) {
      ENZO_FAIL("MFIonizationTestInitialize: either temperature or IEnergy required!\n");
    }
    else {
      RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
      float mp = 1.67262171e-24;    // proton mass [g]
      float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
      float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
      if (RadHydroChemistry == 1) {
	nH = RadHydroDensity*RadHydroHydrogenMassFraction;
	HI = nH*(1.0 - RadHydroInitialFractionHII);
	HII = nH*RadHydroInitialFractionHII;
	ne = HII;
	num_dens = HI + HII + ne;
	mu = RadHydroDensity/num_dens;
      }
      else if (RadHydroChemistry == 3) {
	nH = RadHydroDensity*RadHydroHydrogenMassFraction;
	nHe = RadHydroDensity*(1.0 - RadHydroHydrogenMassFraction);
	HI = nH*(1.0 - RadHydroInitialFractionHII);
	HII = nH*RadHydroInitialFractionHII;
	HeII = nHe*RadHydroInitialFractionHeII;
	HeIII = nHe*RadHydroInitialFractionHeIII;
	HeI = nHe - HeII - HeIII;
	ne = HII + HeII/4.0 + HeIII/2.0;
	num_dens = 0.25*(HeI + HeII + HeIII) + HI + HII + ne;
	mu = RadHydroDensity/num_dens;
      }
      else 
	ENZO_FAIL("MFIonizationTestInitialize: Nchem != {1,3}\n");
      // compute the internal energy
      RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	
    }
  }

  /////////////////
  // Set up the TopGrid as usual
  HierarchyEntry *TempGrid = &TopGrid;
  while (TempGrid != NULL) {
    if (TempGrid->GridData->MFIonizationTestInitializeGrid(
		        RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity,
			RadHydroX1Velocity, RadHydroX2Velocity, RadHydroIEnergy,
			FSRadiation, Radiation1, Radiation2, Radiation3, 
			E1Units, E2Units, E3Units, RadHydroHydrogenMassFraction,
			RadHydroInitialFractionHII, RadHydroInitialFractionHeII,
			RadHydroInitialFractionHeIII, local) == FAIL) 
      ENZO_FAIL("MFIonizationTestInitialize: Grid initialization failure\n");
    TempGrid = TempGrid->NextGridThisLevel;
  }


  /////////////////
  // If requested, refine the grid to the desired level (second go around)
  // NOTE: this only works if the requested patches are nested, i.e. this 
  //       approach can only currently handle one patch per AMR level.
  if (AMRNumberOfInitialPatches > 0  &&  local==1) {

    // Declare, initialize and fill out the LevelArray
    int numlevels = 0;
    for (patch=0; patch<AMRNumberOfInitialPatches; patch++)
      numlevels = max(numlevels,AMRPatchLevel[patch]);
    if (numlevels > MAX_DEPTH_OF_HIERARCHY) 
      ENZO_FAIL("Too many initial patch levels! increase MAX_DEPTH_OF_HIERARCHY!\n");
    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level=0; level<MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    // store the current values for MustRefineRegion, StaticHierarchy
    int MRRMRL_save = MustRefineRegionMinRefinementLevel;
    int CFM_save[MAX_FLAGGING_METHODS];
    for (i=0; i<MAX_FLAGGING_METHODS; i++) {
      CFM_save[i] = CellFlaggingMethod[i];
      CellFlaggingMethod[i] = 0;
    }
    CellFlaggingMethod[0] = 12;
    FLOAT MRRLE_save[MAX_DIMENSION], MRRRE_save[MAX_DIMENSION];
    for (i=0; i<MAX_DIMENSION; i++) {
      MRRLE_save[i] = MustRefineRegionLeftEdge[i];
      MRRRE_save[i] = MustRefineRegionRightEdge[i];
    }
    int SH_save = MetaData.StaticHierarchy;
    MetaData.StaticHierarchy = FALSE;

    // For each patch, adjust the MustRefineRegion parameters and re-call
    // RebuildHierarchy to refine that region as needed.  Then re-initialize
    // the grids after they are.
    for (patch=0; patch<AMRNumberOfInitialPatches; patch++) {

      // set level for this patch
      level = AMRPatchLevel[patch]-1;

      // set MustRefineRegion and level for this patch
      MustRefineRegionMinRefinementLevel = level+1;
      for (i=0; i<MAX_DIMENSION; i++) {
	MustRefineRegionLeftEdge[i]  = AMRPatchLeftEdge[patch][i];
	MustRefineRegionRightEdge[i] = AMRPatchRightEdge[patch][i];
      }
      
      // call RebuildHierarchy and initialize grid on this level
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL)
	ENZO_FAIL("Error in RebuildHierarchy.\n");
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->MFIonizationTestInitializeGrid(
		        RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity,
			RadHydroX1Velocity, RadHydroX2Velocity, RadHydroIEnergy,
			FSRadiation, Radiation1, Radiation2, Radiation3, 
			E1Units, E2Units, E3Units, RadHydroHydrogenMassFraction,
			RadHydroInitialFractionHII, RadHydroInitialFractionHeII,
			RadHydroInitialFractionHeIII, 1) == FAIL) 
	  ENZO_FAIL("MFIonizationTestInitialize: Grid initialization failure\n");
	Temp = Temp->NextGridThisLevel;
      }
    }

    // Loop back from the bottom, restoring the consistency among levels
    for (level=numlevels; level>0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");
	Temp = Temp->NextGridThisLevel;
      }
    }

    // re-set the original values for MustRefineRegion*, StaticHierarchy
    MustRefineRegionMinRefinementLevel = MRRMRL_save;
    for (i=0; i<MAX_FLAGGING_METHODS; i++)  CellFlaggingMethod[i] = CFM_save[i];
    for (i=0; i<MAX_DIMENSION; i++) {
      MustRefineRegionLeftEdge[i] = MRRLE_save[i];
      MustRefineRegionRightEdge[i] = MRRRE_save[i];
    }
    MetaData.StaticHierarchy = SH_save;

  }  // end if (AMRNumberOfInitialPatches > 0  &&  local==0) {


  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = FSRadName;
  DataLabel[BaryonField++] = Rad1Name;
  DataLabel[BaryonField++] = Rad2Name;
  DataLabel[BaryonField++] = Rad3Name;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;
  DataLabel[BaryonField++] = HeIName;
  DataLabel[BaryonField++] = HeIIName;
  DataLabel[BaryonField++] = HeIIIName;

  // if using external chemistry/cooling, set rate labels
  if (RadiativeCooling) {
    DataLabel[BaryonField++] = kphHIName;
    DataLabel[BaryonField++] = gammaName;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      DataLabel[BaryonField++] = kphHeIName;
      DataLabel[BaryonField++] = kphHeIIName;
    }
    if (MultiSpecies > 1)
      DataLabel[BaryonField++] = kdissH2IName;
  }

  // if using the AMRFLDSplit solver, set a field for the emissivity
  if (ImplicitProblem == 6) 
    DataLabel[BaryonField++] = EtaName;


  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  ENZO_FAIL("MFIonizationTestInitialize: TRANSFER must be enabled for this test!");
 
#endif

}
