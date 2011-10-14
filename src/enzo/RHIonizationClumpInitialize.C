/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- CLUMP IONIZATION TEST
/
/  written by: Daniel Reynolds
/  date:       December 2007
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
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);


int RHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Entering RHIonizationClumpInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *DeName    = "Electron_Density";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient gas velocity - free parameter
  //  3. ambient gas temperature
  //  4. ambient radiation energy
  //  5. Hydrogen mass fraction 
  //  6. initial fraction HII
  //  7. Number of chemical species
  //  8. mesh spacing

  const char config_rad_hydro_ionization_clump_defaults[] = 
  "### RAD HYDRO IONIZATION CLUMP DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    RadHydro: {\n"
  "        Velocity = [0.0, 0.0, 0.0];\n"
  "        NumDensityIn         = 0.04;\n"
  "        NumDensityOut        = 0.0002;\n"
  "        TemperatureIn        = 40.0;\n"
  "        TemperatureOut       = 8000.0;\n"
  "        RadiationEnergy      = 1.0e-20;\n"
  "        HydrogenMassFraction = 1.0;\n"
  "        InitialFractionHII   = 0.0;\n"
  "        ClumpCenter = [1.54285e22,1.018281e22,1.018281e22];\n" //cm (5, 3.3,3.3 kpc)
  "        ClumpRadius          = 2.46856e21;\n"  // cm (0.8 kpc)
  "    };\n"
  "    gFLD: {\n"
  "        Chemistry    = 1;\n"
  "        Model        = 1;\n"
  "    };\n"
  "};\n";

  float RadHydroVelocity[MAX_DIMENSION];
  float RadHydroX0Velocity;
  float RadHydroX1Velocity;
  float RadHydroX2Velocity;
  float RadHydroNumDensityIn;
  float RadHydroNumDensityOut;
  float RadHydroTemperatureIn;
  float RadHydroTemperatureOut;
  float RadHydroRadiationEnergy;
  float RadHydroHydrogenMassFraction;
  float RadHydroInitialFractionHII;
  int   RadHydroChemistry;
  int   RadHydroModel;
  float ClumpCenter[MAX_DIMENSION];
  float ClumpRadius;  // cm (0.8 kpc)

  Param.Update(config_rad_hydro_ionization_clump_defaults);

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      // read relevant problem parameters
      Param.GetArray(RadHydroVelocity, "Problem.RadHydro.Velocity");
      Param.GetScalar(RadHydroChemistry, "Problem.gFLD.Chemistry");
      Param.GetScalar(RadHydroModel, "Problem.gFLD.Model");
      Param.GetScalar(RadHydroNumDensityIn, "Problem.RadHydro.NumDensityIn");
      Param.GetScalar(RadHydroNumDensityOut, "Problem.RadHydro.NumDensityOut");
      Param.GetScalar(RadHydroTemperatureIn, "Problem.RadHydro.TemperatureIn");
      Param.GetScalar(RadHydroTemperatureOut, "Problem.RadHydro.TemperatureOut");
      Param.GetScalar(RadHydroRadiationEnergy, "Problem.RadHydro.RadiationEnergy");
      Param.GetScalar(RadHydroInitialFractionHII, "Problem.RadHydro.InitialFractionHII");
      Param.GetArray(ClumpCenter, "Problem.RadHydro.ClumpCenter");
      Param.GetScalar(ClumpRadius, "Problem.RadHydro.ClumpRadius");
      fclose(RHfptr);
    }
  }
  RadHydroX0Velocity = RadHydroVelocity[0];
  RadHydroX1Velocity = RadHydroVelocity[1];
  RadHydroX2Velocity = RadHydroVelocity[2];
  ClumpCenterX = RadHydroCenter[0]
  ClumpCenterY = RadHydroCenter[1]
  ClumpCenterZ = RadHydroCenter[2]

  // ensure that we're performing only Hydrogen chemistry
  if (RadHydroChemistry != 1) 
    ENZO_FAIL("RHIonizationClumpInitialize error: RadHydroChemistry must equal 1!");

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // if temperature specified and not internal energy, perform conversion here
  RadHydroTemperatureIn = max(RadHydroTemperatureIn,MIN_TEMP); // enforce minimum
  RadHydroTemperatureOut = max(RadHydroTemperatureOut,MIN_TEMP); // enforce minimum
  float mp = 1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float nH, HI, HII, ne, num_dens, mu;
  HI = 1.0 - RadHydroInitialFractionHII;
  HII = RadHydroInitialFractionHII;
  ne = HII;
  num_dens = HI + HII + ne;
  mu = 1.0/num_dens;
  // correct mu if using a special model
  if ((RadHydroModel == 4) || (RadHydroModel == 5)) 
    mu = DEFAULT_MU;
  // compute the internal energy
  float RadHydroIEnergyIn  = kb*RadHydroTemperatureIn/mu/mp/(Gamma-1.0);
  float RadHydroIEnergyOut = kb*RadHydroTemperatureOut/mu/mp/(Gamma-1.0);
  
  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RHIonizationClumpInitializeGrid(
			RadHydroChemistry, RadHydroNumDensityIn, 
			RadHydroNumDensityOut, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergyIn, RadHydroIEnergyOut, 
			RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, ClumpCenterX, ClumpCenterY, 
			ClumpCenterZ, ClumpRadius, local) == FAIL) {
      fprintf(stderr, "Error in RHIonizationClumpInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

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
