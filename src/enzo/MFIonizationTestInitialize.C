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
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);




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

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

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
      } // end input from parameter file
      fclose(RHfptr);
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

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->MFIonizationTestInitializeGrid(
		        RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity,
			RadHydroX1Velocity, RadHydroX2Velocity, RadHydroIEnergy,
			FSRadiation, Radiation1, Radiation2, Radiation3, 
			E1Units, E2Units, E3Units, RadHydroHydrogenMassFraction,
			RadHydroInitialFractionHII, RadHydroInitialFractionHeII,
			RadHydroInitialFractionHeIII, local) == FAIL) 
      ENZO_FAIL("MFIonizationTestInitialize: Grid initialization failure\n");
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

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  ENZO_FAIL("MFIonizationTestInitialize: TRANSFER must be enabled for this test!");
 
#endif

}
