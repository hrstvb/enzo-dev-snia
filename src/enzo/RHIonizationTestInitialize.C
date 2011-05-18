/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- IONIZATION TEST
/
/  written by: Daniel Reynolds
/  date:       July 2007
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
#include "LevelHierarchy.h"
#include "TopGridData.h"


/* default constants */
#define DEFAULT_MU 0.6           // mean molecular mass
#define MIN_TEMP 1.0             // minimum temperature [K]
#define MAX_INITIAL_PATCHES 100  // max number of parameter file defined subgrids

// function prototypes
int InitializeRateData(FLOAT Time);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);



int RHIonizationTestInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local)
{
#ifdef TRANSFER
//   if (MyProcessorNumber == ROOT_PROCESSOR)
//     fprintf(stdout,"Entering RHIonizationTestInitialize routine\n");

  // if ParallelRootGridIO is disabled, only perform this routine the first pass 
  // through (when local==0).  This enables initialization of an initially refined 
  // mesh (otherwise things do not get set up properly)
  if (!ParallelRootGridIO) {
    if (local == 1)   return SUCCESS;
    else              local = 1;
  }

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
  char *RadName   = "Grey_Radiation_Energy";
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

  /////////////////
  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient gas velocity - free parameter
  //  3. ambient gas temperature
  //  4. ambient radiation energy
  //  5. Hydrogen mass fraction 
  //  6. initial fraction HII
  //  7. initial fraction HeII
  //  8. initial fraction HeIII
  //  9. Number of chemical species
  // 10. mesh spacing
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroDensity              = 10.0;
  float RadHydroTemperature          = 1.0;
  float RadHydroIEnergy              = -1.0;
  float RadHydroRadiationEnergy      = 10.0;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  float RadHydroInitialFractionHeII  = 0.0;
  float RadHydroInitialFractionHeIII = 0.0;
  int   RadHydroChemistry            = 1;
  int   AMRNumberOfInitialPatches = 0;
  int   AMRPatchLevel[MAX_INITIAL_PATCHES];
  FLOAT AMRPatchLeftEdge[MAX_INITIAL_PATCHES][MAX_DIMENSION];
  FLOAT AMRPatchRightEdge[MAX_INITIAL_PATCHES][MAX_DIMENSION];
  int   AMRPatchDims[MAX_INITIAL_PATCHES][MAX_DIMENSION];
  for (i=0; i<MAX_INITIAL_PATCHES; i++) {
    AMRPatchLevel[i] = 0;
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      AMRPatchLeftEdge[i][dim] = 0.0;
      AMRPatchRightEdge[i][dim] = 0.0;
      AMRPatchDims[i][dim] = 0;
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
	ret += sscanf(line, "RadHydroDensity = %"FSYM, 
		      &RadHydroDensity);
	ret += sscanf(line, "RadHydroTemperature = %"FSYM, 
		      &RadHydroTemperature);
	ret += sscanf(line, "RadHydroIEnergy = %"FSYM, 
		      &RadHydroIEnergy);
	ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
		      &RadHydroRadiationEnergy);
	if (RadHydroChemistry > 0) {
	  ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
			&RadHydroInitialFractionHII);
	  ret += sscanf(line, "RadHydroHFraction = %"FSYM, 
			&RadHydroHydrogenMassFraction);
	}
	if ((RadHydroChemistry == 3) || (MultiSpecies == 1)) {
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


  // if both an initially refined mesh and ParallelRootGridIO are enabled, issue 
  // an error message since that functionality is not currently supported
  if (ParallelRootGridIO && AMRNumberOfInitialPatches) 
    ENZO_FAIL("ParallelRootGridIO and initial refined mesh (AMRNumberOfInitialPatches) incompatible!\n");


  // set grid dimensions of each patch
  float dx;
  for (patch=0; patch<AMRNumberOfInitialPatches; patch++) 
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      dx = 1.0 / float(MetaData.TopGridDims[dim]) / POW(RefineBy, AMRPatchLevel[patch]);
      AMRPatchDims[patch][dim] = 
	nint((AMRPatchRightEdge[patch][dim] - AMRPatchLeftEdge[patch][dim]) / 
	     (DomainRightEdge[dim] - DomainLeftEdge[dim]) / dx);		
    }  


  // output requested initial AMR hierarchy structure
  if (debug && local && AMRNumberOfInitialPatches) {
    printf("\n  Initializing with prescribed AMR hierarchy:\n");
    printf("    AMRNumberOfInitialPatches = %"ISYM"\n\n",
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
      printf("    AMRPatchDims[%"ISYM"] = %"ISYM" %"ISYM" %"ISYM"\n\n", patch, 
	     AMRPatchDims[patch][0], AMRPatchDims[patch][1], AMRPatchDims[patch][2]);
    }
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.\n");

  // if temperature specified and not internal energy, perform conversion here
  if (RadHydroIEnergy == -1.0) {
    if (RadHydroTemperature == -1.0) {
      ENZO_FAIL("Initialize error: either temperature or IEnergy required!\n");
    }
    else {
      RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
      float mp = 1.67262171e-24;    // proton mass [g]
      float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
      float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
      if (RadHydroChemistry == 0) 
	mu = DEFAULT_MU;
      else if (RadHydroChemistry == 1) {
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
      else {
	ENZO_FAIL("Initialize error: NChem != {0,1,3}\n");
      }
      // compute the internal energy
      RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	
    }
  }

  /////////////////
  // Set up the TopGrid as usual
  HierarchyEntry *TempGrid = &TopGrid;
  while (TempGrid != NULL) {
    if (TempGrid->GridData->RHIonizationTestInitializeGrid(
		        RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, RadHydroIEnergy, 
			RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, RadHydroInitialFractionHeII, 
			RadHydroInitialFractionHeIII, local) == FAIL) {
      ENZO_FAIL("Error in RHIonizationTestInitializeGrid!\n");
    }
    TempGrid = TempGrid->NextGridThisLevel;
  }
  

  /////////////////
  // If requested, add patches -- currently requires strict nesting of subgrids, 
  // with one patch per level
  if (AMRNumberOfInitialPatches > 0  &&  local==1) {

    // determine how many levels we'll fill in
    int numlevels = 0;
    for (patch=0; patch<AMRNumberOfInitialPatches; patch++)
      numlevels = max(numlevels,AMRPatchLevel[patch]);
    if (numlevels > MAX_DEPTH_OF_HIERARCHY) 
      ENZO_FAIL("Too many initial patch levels! increase MAX_DEPTH_OF_HIERARCHY!\n");
    if (numlevels != AMRNumberOfInitialPatches)
      ENZO_FAIL("Error: patches must be strictly nested, one per level (for now)!\n");

    if (numlevels > 0) {
      
      // data for setting up initially-refined regions
      HierarchyEntry **Subgrid = new HierarchyEntry*[numlevels];
      for (patch=0; patch<numlevels; patch++) 
	Subgrid[patch] = new HierarchyEntry;

      for (patch=0; patch<numlevels; patch++) {
	
	int level = AMRPatchLevel[patch];
	if (AMRPatchDims[patch][0] > 0) {

	  for (dim = 0; dim < MetaData.TopGridRank; dim++)
	    AMRPatchDims[patch][dim] += 2*DEFAULT_GHOST_ZONES;
	  
	  // Insert into AMR hierarchy
	  if (level == 1) {
	    Subgrid[patch]->NextGridThisLevel = TopGrid.NextGridNextLevel;
	    TopGrid.NextGridNextLevel = Subgrid[patch];
	    Subgrid[patch]->ParentGrid = &TopGrid;
	  } else {
	    Subgrid[patch]->NextGridThisLevel = NULL;
	    Subgrid[patch]->ParentGrid = Subgrid[patch-1];
	  }
	  if (level == numlevels) 
	    Subgrid[patch]->NextGridNextLevel = NULL;
	  else 
	    Subgrid[patch]->NextGridNextLevel = Subgrid[patch+1];
	  
	  // Create grid
	  Subgrid[patch]->GridData = new grid;
	  Subgrid[patch]->GridData->InheritProperties(Subgrid[patch]->ParentGrid->GridData);
	  Subgrid[patch]->GridData->PrepareGrid(MetaData.TopGridRank, 
						AMRPatchDims[patch],
						AMRPatchLeftEdge[patch],
						AMRPatchRightEdge[patch], 0);

	  if (Subgrid[patch]->GridData->RHIonizationTestInitializeGrid(
		        RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, RadHydroIEnergy, 
			RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, RadHydroInitialFractionHeII, 
			RadHydroInitialFractionHeIII, local) == FAIL)
		ENZO_FAIL("Error in RHIonizationTestInitializeGrid.\n");
	} // endif (AMRPatchDims[patch][0] > 0)
      } // endfor (patch)
    } // endif (numlevels > 0)


    // All patches and data have been set up fully at this time.  However, we'll now 
    // build a temporary LevelArray so that we can ensure "consistency" between levels
    
    // Declare, initialize and fill out the LevelArray
    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level=0; level<MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    // Loop back from the bottom, restoring the consistency among levels
    for (level=MaximumRefinementLevel; level>0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
			      *LevelArray[level-1]->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
	}
	//	Temp->GridData->ClearBoundaryFluxes();
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end if (AMRNumberOfInitialPatches > 0  &&  local==1) {


//   // output initial grid properties
//   printf("\nRHIonizationTestInitialize: grids and their properties:\n");
//   for (HierarchyEntry* Temp=&TopGrid; Temp; Temp=Temp->NextGridNextLevel) {
//     for (HierarchyEntry* Temp2=Temp; Temp; Temp2=Temp2->NextGridThisLevel) {
//       if (Temp2 == NULL) break;
//       printf("   p%"ISYM": Temp = %p, Temp2 = %p\n",
// 	     MyProcessorNumber,Temp,Temp2);
//       printf("   p%"ISYM": grid %p: NumberOfBaryonFields = %"ISYM", location = (%g,%g) x (%g,%g) x (%g,%g)\n",
// 	     MyProcessorNumber, Temp2, 
// 	     Temp2->GridData->ReturnNumberOfBaryonFields(), 
// 	     Temp2->GridData->GetGridLeftEdge(0), 
// 	     Temp2->GridData->GetGridRightEdge(0), 
// 	     Temp2->GridData->GetGridLeftEdge(1), 
// 	     Temp2->GridData->GetGridRightEdge(1), 
// 	     Temp2->GridData->GetGridLeftEdge(2), 
// 	     Temp2->GridData->GetGridRightEdge(2));
//       printf("   p%"ISYM": grid %p: NGNL = %p, NGTL = %p, PG = %p\n\n",
// 	     MyProcessorNumber,Temp2, 
// 	     Temp2->NextGridNextLevel, Temp2->NextGridThisLevel, Temp2->ParentGrid);
//     }
//   }  // end loop over grids


  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName;
  if (RadHydroChemistry > 0) {
    DataLabel[BaryonField++] = DeName;
    DataLabel[BaryonField++] = HIName;
    DataLabel[BaryonField++] = HIIName;
  }
  if ((RadHydroChemistry == 3) || (MultiSpecies > 0)) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

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

  ENZO_FAIL("Error: TRANSFER must be enabled for this test!\n");
 
#endif

}
