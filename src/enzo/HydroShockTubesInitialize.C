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
#include "LevelHierarchy.h"
#include "TopGridData.h"

/* Set default parameter values. */

const char config_hydro_shocktubes_defaults[] =
"### HYDRO SHOCKTUBES DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    HydroShocktubes: {\n"
"        RefineAtStart		= FALSE;\n"
"        InitialDiscontinuity 	= 0.5;\n"
"        SecondDiscontinuity 	= 0.5;\n"
"\n"
"        LeftDensity 		= 1.0;\n"
"	 RightDensity 		= 1.0;\n"
"        CenterDensity 		= 1.0;\n"
"\n"
"        LeftPressure 		= 1.0;\n" 
"        RightPressure 		= 1.0;\n"
"        CenterPressure 	= 1.0;\n"
"\n"
"        LeftVelocity   = [0.0,0.0,0.0];\n"
"        RightVelocity  = [0.0,0.0,0.0];\n"
"        CenterVelocity = [1.0,1.0,1.0];\n"
"    };\n"
"};\n";



int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int HydroShockTubesInitialize(FILE *fptr, FILE *Outfptr, 
			      HierarchyEntry &TopGrid, TopGridData &MetaData) 
{
  char *DensName = "Density";
  char *PresName = "Pressure";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  float  InitialDiscontinuity, SecondDiscontinuity,
    LeftDensity, RightDensity, CenterDensity, 
    LeftPressure, RightPressure, CenterPressure;
  float LeftVelocity[3];
  float CenterVelocity[3];
  float RightVelocity[3];
 
  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_hydro_shocktubes_defaults);
 

  /* read parameters */
  
  Param.GetScalar(HydroShockTubesRefineAtStart,
		"Problem.HydroShocktube.RefineAtStart");
  Param.GetScalar(HydroShockTubesInitialDiscontinuity,
		"Problem.HydroShocktube.InitialDiscontinuity");
  Param.GetScalar(HydroShockTubesSecondDiscontinuity,
		"Problem.HydroShocktube.SecondDiscontinuity");

  Param.GetArray(HydroShockTubesLeftVelocity,
		"Problem.HydroShocktube.LeftVelocity");
  Param.GetScalar(HydroShockTubesLeftPressure,
		"Problem.HydroShocktube.LeftPressure");
  Param.GetScalar(HydroShockTubesLeftDensity,
		"Problem.HydroShocktube.LeftDensity");

  Param.GetArray(HydroShockTubesRightVelocity, 
		"Probem.HydroShocktube.RightVelocity");
  Param.GetScalar(HydroShockTubesRightPressure,
		"Probem.HydroShocktube.RightPressure");
  Param.GetScalar(HydroShockTubesRightDensity,
		"Probem.HydroShocktube.RightDensity");

  Param.GetArray(HydroShockTubesCenterVelocity,
		"Probem.HydroShocktube.CenterVelocity");
  Param.GetScalar(HydroShockTubesCenterPressure, 
		"Probem.HydroShocktube.CenterPressure");
  Param.GetScalar(HydroShockTubesCenterDensity,
                "Probem.HydroShocktube.CenterDensity");

  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set up grid */

  TopGrid.GridData->
    HydroShockTubesInitializeGrid(InitialDiscontinuity, 
				  SecondDiscontinuity,
				  LeftDensity, RightDensity, CenterDensity,
				  LeftVelocity[0],  RightVelocity[0], CenterVelocity[0],
				  LeftVelocity[1],  RightVelocity[1], CenterVelocity[1],
				  LeftVelocity[2],  RightVelocity[2], CenterVelocity[2],
				  LeftPressure,   RightPressure,  CenterPressure);

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (RefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      printf("In level %"ISYM"\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	Temp->GridData->
	  HydroShockTubesInitializeGrid
	  (InitialDiscontinuity, SecondDiscontinuity,
	   LeftDensity, RightDensity, CenterDensity,
	   LeftVelocity[0],  RightVelocity[0], CenterVelocity[0],
	   LeftVelocity[1],  RightVelocity[1], CenterVelocity[1],
	   LeftVelocity[2],  RightVelocity[2], CenterVelocity[2],
	   LeftPressure,   RightPressure,  CenterPressure);
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels


    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

    //WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
    //       &TopGrid, MetaData, Exterior, -1);

  } // end: if (RefineAtStart)


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "HydroShockTubesRefineAtStart        = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "HydroShockTubesInitialDiscontinuity = %"FSYM"\n",
	    InitialDiscontinuity);
    fprintf(Outfptr, "HydroShockTubesLeftDensity          = %"FSYM"\n",
	    LeftDensity);
    fprintf(Outfptr, "HydroShockTubesRightDensity         = %"FSYM"\n",
	    RightDensity);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityX        = %"FSYM"\n",
	    LeftVelocity[0]);
    fprintf(Outfptr, "HydroShockTubesRightVelocityX       = %"FSYM"\n",
            RightVelocity[0]);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityY        = %"FSYM"\n",
	    LeftVelocity[1]);
    fprintf(Outfptr, "HydroShockTubesRightVelocityY       = %"FSYM"\n",
            RightVelocity[1]);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityZ        = %"FSYM"\n",
	    LeftVelocity[2]);
    fprintf(Outfptr, "HydroShockTubesRightVelocityZ       = %"FSYM"\n",
            RightVelocity[2]);
    fprintf(Outfptr, "HydroShockTubesLeftPressure         = %"FSYM"\n",
            LeftPressure);
    fprintf(Outfptr, "HydroShockTubesRightPressure        = %"FSYM"\n",
            RightPressure);

    fprintf(Outfptr, "HydroShockTubesSecondDiscontinuity = %"FSYM"\n",
	    SecondDiscontinuity);
    fprintf(Outfptr, "HydroShockTubesCenterDensity       = %"FSYM"\n",
	    CenterDensity);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityX     = %"FSYM"\n",
	    CenterVelocity[0]);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityY     = %"FSYM"\n",
	    CenterVelocity[1]);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityZ     = %"FSYM"\n",
	    CenterVelocity[2]);
    fprintf(Outfptr, "HydroShockTubesCenterPressure      = %"FSYM"\n",
	    CenterPressure);

  }
  //return FAIL;
  return SUCCESS;

}

