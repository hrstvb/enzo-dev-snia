/***********************************************************************
/
/  INITIALIZE A ZELDOVICH PANCAKE
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//

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
#include "CosmologyParameters.h"



 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int ZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";


  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int ret;
 
  /* Error check. */
 
  if (!ComovingCoordinates) {
    ENZO_FAIL("ComovingCoordinates must be TRUE!\n");
  }
 
  if (!SelfGravity)
    fprintf(stderr, "ZeldovichPancake: gravity is off!?!\n");
  if (CellFlaggingMethod[0] < 2)
    fprintf(stderr, "ZeldovichPancake: check CellFlaggingMethod.\n");
 
  /* set default parameters */
  const char config_zeldovich_pancake_defaults[] = 
  "### ZELDOVICH PANCAKE DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    ZeldovichPancake: {\n"
  "        Direction = 0;\n"
  "        CentralOffset      = 0.0;\n"
  "        OmegaBaryonNow     = 1.0;\n"
  "        OmegaCDMNow        = 0.0;\n"
  "        CollapseRedshift   = 1.0;\n"
  "        InitialTemperature = 100;\n"
  "        InitialUniformBField =  [0.0, 0.0, 0.0];\n"
  "    };\n"
  "};\n"; 
  int   ZeldovichPancakeDirection;    // along the x-axis
  float ZeldovichPancakeCentralOffset;  // no offset
  float ZeldovichPancakeOmegaBaryonNow;  // standard
  float ZeldovichPancakeOmegaCDMNow;  // no dark matter
  float ZeldovichPancakeCollapseRedshift;  // free parameter
  float ZeldovichPancakeInitialTemperature;  // whatever
  float ZeldovichPancakeInitialUniformBField[MAX_DIMENSION];  // in Gauss

  Param.Update(config_zeldovich_pancake_defaults);
 
    /* read parameters */
 
  Param.GetScalar(ZeldovichPancakeDirection, "Problem.ZeldovichPancake.Direction");
  Param.GetScalar(ZeldovichPancakeCentralOffset, "Problem.ZeldovichPancake.CentralOffset");
  Param.GetScalar(ZeldovichPancakeOmegaBaryonNow, "Problem.ZeldovichPancake.OmegaBaryonNow");
  Param.GetScalar(ZeldovichPancakeOmegaCDMNow, "Problem.ZeldovichPancake.OmegaCDMNow");
  Param.GetScalar(ZeldovichPancakeCollapseRedshift, "Problem.ZeldovichPancake.CollapseRedshift");
  Param.GetScalar(ZeldovichPancakeInitialTemperature, "Problem.ZeldovichPancake.InitialTemperature");
  Param.GetArray(ZeldovichPancakeInitialUniformBField, "Problem.ZeldovichPancake.InitialUniformBField");



  /* Convert from Gauss */
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1, PressureUnits=1.,MagneticUnits=1., a=1,dadt=0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  PressureUnits = DensityUnits * (LengthUnits/TimeUnits)*(LengthUnits/TimeUnits);
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);

  for (int dim = 0; dim < MAX_DIMENSION; dim++) 
    ZeldovichPancakeInitialUniformBField[dim] /= MagneticUnits;

 
  /* set up grid */
 
  if (TopGrid.GridData->ZeldovichPancakeInitializeGrid(
					  ZeldovichPancakeDirection,
					  ZeldovichPancakeCentralOffset,
					  ZeldovichPancakeOmegaBaryonNow,
					  ZeldovichPancakeOmegaCDMNow,
					  ZeldovichPancakeCollapseRedshift,
					  ZeldovichPancakeInitialTemperature,
					  ZeldovichPancakeInitialUniformBField
						       ) == FAIL) {
    ENZO_FAIL("Error in ZeldovichPancakeInitializeGrid.\n");
  }
 
  /* set up field names and units */
 
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = Vel1Name;
  if (MetaData.TopGridRank > 1 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    DataLabel[i++] = Vel2Name;
  if (MetaData.TopGridRank > 2 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    DataLabel[i++] = Vel3Name;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
    if(UseDivergenceCleaning){
      DataLabel[i++] = Phi_pName;
      DataLabel[i++] = DebugName;
    }
  }
 
  for (int count = 0; count < i; count++)
    DataUnits[count] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "ZeldovichPancakeDirection          = %"ISYM"\n",
	    ZeldovichPancakeDirection);
    fprintf(Outfptr, "ZeldovichPancakeCentralOffset      = %"FSYM"\n",
	    ZeldovichPancakeCentralOffset);
    fprintf(Outfptr, "ZeldovichPancakeOmegaBaryonNow     = %"FSYM"\n",
	    ZeldovichPancakeOmegaBaryonNow);
    fprintf(Outfptr, "ZeldovichPancakeOmegaCDMNow        = %"FSYM"\n",
	    ZeldovichPancakeOmegaCDMNow);
    fprintf(Outfptr, "ZeldovichPancakeCollapseRedshift   = %"FSYM"\n",
	    ZeldovichPancakeCollapseRedshift);
    fprintf(Outfptr, "ZeldovichPancakeInitialTemperature = %"FSYM"\n\n",
	    ZeldovichPancakeInitialTemperature);

    for (int dim = 0; dim < MAX_DIMENSION; dim++) 
      ZeldovichPancakeInitialUniformBField[dim] *= MagneticUnits;
    fprintf(Outfptr, "ZeldovichPancakeInitialUniformBField = ");
    WriteListOfFloats(Outfptr, 3, ZeldovichPancakeInitialUniformBField);

  }
 
  return SUCCESS;
}
