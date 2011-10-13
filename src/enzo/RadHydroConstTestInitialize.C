/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- EVOLVE A CONSTANT FIELD
/
/  written by: Daniel Reynolds
/  date:       November, 2006
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




int RadHydroConstTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering RadHydroConstTestInitialize routine\n");

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
  //  7. initial fraction HeII
  //  8. initial fraction HeIII
  //  9. Number of chemical species
  // 10. mesh spacing
  float RadHydroX0Velocity;
  float RadHydroX1Velocity;
  float RadHydroX2Velocity;
  float RadHydroDensity;
  float RadHydroTemperature;
  float RadHydroIEnergy;
  float RadHydroRadiationEnergy;
  float RadHydroHydrogenMassFraction;
  float RadHydroInitialFractionHII;
  float RadHydroInitialFractionHeII;
  float RadHydroInitialFractionHeIII;
  int   RadHydroChemistry;
  int   RadHydroModel;

  const char config_rad_hydro_const_test_defaults[] = 
  "### RAD HYDRO CONST TEST DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    RadHydro: {\n"
  "        Velocity             = [0.0, 0.0, 0.0];\n"
  "        Density              = 10.0;\n"
  "        Temperature          = 1.0;\n"
  "        IEnergy              = -1.0;\n"
  "        RadiationEnergy      = 10.0;\n"
  "        HydrogenMassFraction = 1.0;\n"
  "        InitialFractionHII   = 0.0;\n"
  "        InitialFractionHeII  = 0.0;\n"
  "        InitialFractionHeIII = 0.0;\n"
  "        Chemistry            = 1;\n"
  "        Model                = 1;\n"
  "    };\n"
  "};\n";

  Param.Update(config_rad_hydro_const_test_defaults);

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      // read relevant problem parameters
      Param.GetArray(RadHydroVelocity, "Problem.RadHydro.Velocity");
      Param.GetScalar(RadHydroChemistry, "Problem.RadHydro.Chemistry");
      Param.GetScalar(RadHydroModel, "Problem.RadHydro.Model");
      Param.GetScalar(RadHydroDensity, "Problem.RadHydro.Density");
      Param.GetScalar(RadHydroTemperature, "Problem.RadHydro.Temperature");
      Param.GetScalar(RadHydroIEnergy, "Problem.RadHydro.IEnergy");
      Param.GetScalar(RadHydroRadiationEnergy, "Problem.RadHydro.RadiationEnergy");
      if (RadHydroChemistry > 0)
	Param.GetScalar(RadHydroInitialFractionHII, "Problem.RadHydro.InitialFractionHII");
      if (RadHydroChemistry > 1) {
	Param.GetScalar(RadHydroHydrogenMassFraction, "Problem.RadHydro.HFraction");
	Param.GetScalar(RadHydroInitialFractionHeII, "Problem.RadHydro.InitialFractionHeII");
	Param.GetScalar(RadHydroInitialFractionHeIII, "Problem.RadHydro.InitialFractionHeIII");
      }
      fclose(RHfptr);
    }
  }
  RadHydroX0Velocity = RadHydroVelocity[0];
  RadHydroX1Velocity = RadHydroVelocity[1];
  RadHydroX2Velocity = RadHydroVelocity[2];



  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // if temperature specified and not internal energy, perform conversion here
  if (RadHydroIEnergy == -1.0) {
    if (RadHydroTemperature == -1.0) {
      fprintf(stderr,"Initialize error: either temperature or IEnergy required!\n");
      return FAIL;
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
	fprintf(stderr,"Initialize error: NChem != {0,1,3}\n");
	return FAIL;	
      }
      // correct mu if using a special model
      if ((RadHydroModel == 4) || (RadHydroModel == 5)) 
	mu = DEFAULT_MU;
      // compute the internal energy
      RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	
    }
  }

  // set up the grid(s) on this level
  printf("RadHydroConstTestInitialize: calling grid initializer\n");
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RadHydroConstTestInitializeGrid(RadHydroChemistry, 
			RadHydroDensity, RadHydroX0Velocity, RadHydroX1Velocity, 
			RadHydroX2Velocity, RadHydroIEnergy, 
			RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
                        RadHydroInitialFractionHII, RadHydroInitialFractionHeII, 
                        RadHydroInitialFractionHeIII, local) == FAIL) {
      fprintf(stderr, "Error in RadHydroConstTestInitializeGrid.\n");
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
  if (RadHydroChemistry > 0) {
    DataLabel[BaryonField++] = DeName;
    DataLabel[BaryonField++] = HIName;
    DataLabel[BaryonField++] = HIIName;
  }
  if (RadHydroChemistry == 3) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }
  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}
