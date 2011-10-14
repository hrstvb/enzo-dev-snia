/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Daniel Reynolds
/  date:       May 2008
/  modified1:
/
/  PURPOSE:
/    Set up a uniform cosmological HII region ionization test
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
 
// default constants
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]
#define MAX_INITIAL_GRIDS 10


// function prototypes
int InitializeRateData(FLOAT Time);

 
int CosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, 
			      TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Entering CosmoIonizationInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *GEName    = "GasEnergy";
  char *Vel1Name  = "x-velocity";
  char *Vel2Name  = "y-velocity";
  char *Vel3Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *DeName    = "Electron_Density";
 
  // local declarations
  char  line[MAX_LINE_LENGTH];
  int   dim, ret;

  const char config_cosmo_ionization_defaults[] = 
  "### COSMO IONIZATION CLUMP DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    RadHydro: {\n"
  "        Velocity = [0.0, 0.0, 0.0];\n"
  "        Temperature          = 10000.0;\n"       // [K]
  "        RadiationEnergy      = 1.0e-32;\n"
  "        InitialFractionHII   = 0.0;\n"
  "        OmegaBaryonNow       = 0.2;\n"
  "    };\n"
  "    gFLD: {\n"
  "        Chemistry    = 1;\n"
  "        Model        = 1;\n"
  "    };\n"
  "};\n";
 
  // Setup and parameters:
  float RadHydroVelocity[MAX_DIMENSION];
  float RadHydroX0Velocity;
  float RadHydroX1Velocity;
  float RadHydroX2Velocity;
  float RadHydroTemperature;       // [K]
  float RadHydroRadiationEnergy;
  float RadHydroInitialFractionHII;
  float RadHydroOmegaBaryonNow;
  int   RadHydroChemistry;
  int   RadHydroModel;

  Param.Update(config_cosmo_ionization_defaults);

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      // read relevant problem parameters
      Param.GetArray(RadHydroVelocity, "Problem.RadHydro.Velocity");
      Param.GetScalar(RadHydroChemistry, "Problem.gFLD.Chemistry");
      Param.GetScalar(RadHydroModel, "Problem.gFLD.Model");
      Param.GetScalar(RadHydroTemperature, "Problem.RadHydro.Temperature");
      Param.GetScalar(RadHydroRadiationEnergy, "Problem.RadHydro.RadiationEnergy");
      Param.GetScalar(RadHydroInitialFractionHII, "Problem.RadHydro.InitialFractionHII");
      Param.GetScalar(RadHydroOmegaBaryonNow, "Problem.RadHydro.OmegaBaryonNow");
      fclose(RHfptr);
    }
  }
  RadHydroX0Velocity = RadHydroVelocity[0];
  RadHydroX1Velocity = RadHydroVelocity[1];
  RadHydroX2Velocity = RadHydroVelocity[2];

  // require Hydrogen only chemistry
  if (RadHydroChemistry != 1) 
    ENZO_FAIL("CosmoIonizationInitialize error: RadHydroChemistry must equal 1!");

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // convert input temperature to internal energy
  RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
  float mp  =1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float HI = 1.0 - RadHydroInitialFractionHII;
  float HII = RadHydroInitialFractionHII;
  float ne = HII;
  float num_dens = HI + HII + ne;
  float mu = 1.0/num_dens;
  // correct mu if using a special model
  if ((RadHydroModel == 4) || (RadHydroModel == 5)) 
    mu = DEFAULT_MU;
  // compute the internal energy
  float RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->CosmoIonizationInitializeGrid(
		        RadHydroChemistry, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergy, RadHydroRadiationEnergy, 
			RadHydroInitialFractionHII, 
			RadHydroOmegaBaryonNow, local) == FAIL) {
      fprintf(stderr, "Error in CosmoIonizationInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }


  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = GEName;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = Vel3Name;
  DataLabel[BaryonField++] = RadName;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}

