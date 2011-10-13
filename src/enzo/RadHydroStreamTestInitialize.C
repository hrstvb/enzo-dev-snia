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


// function prototypes
int InitializeRateData(FLOAT Time);




int RadHydroStreamTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				 TopGridData &MetaData, int local)
{

#ifdef TRANSFER

  if (debug)
    fprintf(stdout,"Entering RadHydroStreamTestInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";

  // local declarations
  int dim;

  printf("Setting up problem with rank = %"ISYM"\n",MetaData.TopGridRank);

  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient radiation energy
  //  3. initial time step size
  //  4. dimension for streaming test problem {0,1,2}
  //  5. direction for streaming radiation {0:l->r, 1:r->l}

  const char config_rad_hydro_stream_test_defaults[] = 
  "### RAD HYDRO STREAM DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    RadHydro: {\n"
  "        Density   = 1.0;\n"
  "        RadEnergy = 1.0e-10;\n"
  "        RadStreamDim        = 0;\n"
  "        RadStreamDir        = 0;\n"
  "    };\n"
  "};\n"; 

  float RadHydroDensity;
  float RadHydroRadEnergy;
  int RadStreamDim;
  int RadStreamDir;

  Param.Update(config_rad_hydro_stream_test_defaults);

  // overwrite parameters from RadHydroParamFile file, if it exists
  char line[MAX_LINE_LENGTH];
  int  ret;
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      // read relevant problem parameters
      Param.GetScalar(RadHydroDensity, "Problem.RadHydro.Density");
      Param.GetScalar(RadHydroRadEnergy, "Problem.RadHydro.RadEnergy");
      Param.GetScalar(RadStreamDim, "Problem.RadHydro.RadStreamDim");
      Param.GetScalar(RadStreamDir, "Problem.RadHydro.RadStreamDir");
      fclose(RHfptr);
    }
  }

  // ensure that streaming dimension is active for this rank
  if (RadStreamDim >= MetaData.TopGridRank) {
    fprintf(stderr,"RadStreamDim = %"ISYM" illegal for rank = %"ISYM"!\n",
	    RadStreamDim,MetaData.TopGridRank);
    return FAIL;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }


  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->RadHydroStreamTestInitializeGrid(RadHydroDensity, 
		        RadHydroRadEnergy, RadStreamDim, RadStreamDir, 
			local) == FAIL) {
      fprintf(stderr, "Error in RadHydroStreamTestInitializeGrid.\n");
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

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#endif

}
