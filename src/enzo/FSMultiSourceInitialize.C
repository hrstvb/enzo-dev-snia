/***********************************************************************
/
/  INITIALIZE FREE-STREAMING RADIATION TEST -- INITIALLY ZERO FIELD 
/  WITH MULTIPLE RANDOMLY-LOCATED SOURCES IN EACH SUBDOMAIN
/
/  written by: Daniel Reynolds
/  date:       June 2009
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

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

const char config_cosmology_simulation_defaults[] = 
"### FS MULTI SOURCE DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    FSMultiSource: {\n"
"         Density       = 10.0;\n"
"         TEnergy       = 1.0;\n"
"         RadiationEnergy  =10.0;\n"
"         Velocity  = [0.0,0.0,0.0];\n"
"    };\n"
"};\n";

/* default constants */
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);




int FSMultiSourceInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			    TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stdout,"Entering FSMultiSourceInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "FS_Radiation";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  float Density              = 10.0;
  float X0Velocity           = 0.0;
  float X1Velocity           = 0.0;
  float X2Velocity           = 0.0;
  float TEnergy              = 1.0;
  float RadiationEnergy      = 10.0;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
  float temparray[3];
  Param.GetArray(temparray,"Problem.FSMultiSource.Velocity");
  X0Velocity=temparray[0];
  X1Velocity=temparray[1];
  X2Velocity=temparray[2];

	Param.GetScalar(FSMultiSourceDensity,"Problem.FSMultiSource.Density");
  Param.GetScalar(FSMultiSourceTEnergy,"Problem.FSMultiSource.TEnergy");
  Param.GetScalar(FSMultiSourceRadiationEnergy,"Problem.FSMultiSource.RadiationEnergy");
      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  // set up the grid(s) on this level
  if (debug)
    printf("FSMultiSourceInitialize: calling grid initializer\n");
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->FSMultiSourceInitializeGrid(Density, X0Velocity, 
			X1Velocity, X2Velocity, TEnergy, 
			RadiationEnergy, local) == FAIL) {
      fprintf(stderr, "Error in Grid_FSMultiSourceInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}
