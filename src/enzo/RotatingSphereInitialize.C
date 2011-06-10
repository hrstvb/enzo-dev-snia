/***********************************************************************
/
/  Rotating Sphere Test Problem
/
/  written by: Brian O'Shea
/  date:       May 2011
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the enzo parameter file
 
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
 
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int RotatingSphereInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
  if(debug){
    printf("Entering RotatingSphereInitialize\n");
    fflush(stdout);
  }

  /* initialize field names */

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *MetalName = "Metal_Density";

/* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, level;
 
  /* make sure it is 3D - this particular problem type breaks down for dims < 3 */
 
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do RotatingSphere in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* position information */
 
  FLOAT RotatingSphereSubgridLeft[MAX_DIMENSION], RotatingSphereSubgridRight[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT RotatingSphereCenterPosition[MAX_DIMENSION];

  /* some declarations and default settings for problem-specific parameters */
 
  for(i=0; i<MAX_DIMENSION; i++)
    RotatingSphereCenterPosition[i] = 0.5;  // right in the middle of the box

  float RotatingSphereVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float RotatingSphereBField[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  FLOAT RotatingSphereCoreRadius = 0.1;
  float RotatingSphereLambda = 0.05;
  float RotatingSphereCentralDensity = 100.0;
  float RotatingSphereCentralTemperature = 1000.0;
  float RotatingSphereExteriorDensity = 1.0;
  float RotatingSphereTotalEnergy = 1.0;
  float Pi                      = 3.14159;
  int RotatingSphereRefineAtStart = FALSE;

  /* set no subgrids by default. */
 
  RotatingSphereSubgridLeft[0] = RotatingSphereSubgridLeft[1] = 
    RotatingSphereSubgridLeft[2] = 0.0;    // start of subgrid(s)

  RotatingSphereSubgridRight[0] = RotatingSphereSubgridRight[1] = 
    RotatingSphereSubgridRight[2] = 0.0;    // end of subgrid(s)

  /* read input from enzo parameter file - just the stuff for this problem type and
     any TestProblem data. */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters specifically for radiating shock problem*/

    ret += sscanf(line, "RotatingSphereCentralDensity  = %"FSYM, &RotatingSphereCentralDensity);
    ret += sscanf(line, "RotatingSphereCentralTemperature  = %"FSYM, &RotatingSphereCentralTemperature);
    ret += sscanf(line, "RotatingSphereSubgridLeft = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereSubgridLeft,RotatingSphereSubgridLeft+1,RotatingSphereSubgridLeft+2);
    ret += sscanf(line, "RotatingSphereSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereSubgridRight,RotatingSphereSubgridRight+1,RotatingSphereSubgridRight+2);
    ret += sscanf(line, "RotatingSphereLambda = %"FSYM,
		        &RotatingSphereLambda);

    ret += sscanf(line, "RotatingSphereTotalEnergy = %"FSYM,
		        &RotatingSphereTotalEnergy);

    ret += sscanf(line, "RotatingSphereRefineAtStart = %"ISYM,
		        &RotatingSphereRefineAtStart);

    ret += sscanf(line, "RotatingSphereCoreRadius = %"PSYM,
		        &RotatingSphereCoreRadius);
    ret += sscanf(line, "RotatingSphereCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereCenterPosition, RotatingSphereCenterPosition+1,
		  RotatingSphereCenterPosition+2);

    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "TestProblemDeuteriumToHydrogenRatio = %"FSYM, &TestProblemData.DeuteriumToHydrogenRatio);

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemInitialHMFraction  = %"FSYM, &TestProblemData.HM_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IFraction  = %"FSYM, &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IIFraction  = %"FSYM, &TestProblemData.H2II_Fraction);
    ret += sscanf(line, "TestProblemInitialDIFraction  = %"FSYM, &TestProblemData.DI_Fraction);
    ret += sscanf(line, "TestProblemInitialDIIFraction  = %"FSYM, &TestProblemData.DII_Fraction);
    ret += sscanf(line, "TestProblemInitialHDIFraction  = %"FSYM, &TestProblemData.HDI_Fraction);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && (strstr(line, "RotatingSphere") || strstr(line, "TestProblem")) &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	 "*** warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file
 
  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

  /* ----------------------------------------------------------------------- */


  if (TopGrid.GridData->InitializeUniformGrid(RotatingSphereExteriorDensity,
					      RotatingSphereTotalEnergy,
					      RotatingSphereTotalEnergy,
					      RotatingSphereVelocity,
					      RotatingSphereBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  if (TopGrid.GridData->RotatingSphereInitializeGrid(RotatingSphereCoreRadius,
						     RotatingSphereCenterPosition,
						     RotatingSphereLambda,
						     RotatingSphereCentralDensity,
						     RotatingSphereCentralTemperature,
						     RotatingSphereExteriorDensity) == FAIL) {
    ENZO_FAIL("Error in RotatingSphereInitializeGrid.");
  }
 
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *= (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  printf("%e %e\n",MinimumMassForRefinement[0],MinimumOverDensityForRefinement[0]);

  /* If requested, refine the grid to the desired level. */

  if (RotatingSphereRefineAtStart) {

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

      printf("Going to create a new grid!\n");
      fflush(stdout);

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	Temp->GridData->RotatingSphereInitializeGrid(RotatingSphereCoreRadius,
						     RotatingSphereCenterPosition,
						     RotatingSphereLambda,
						     RotatingSphereCentralDensity,
						     RotatingSphereCentralTemperature,
						     RotatingSphereExteriorDensity) ;
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

  } // end: if (RotatingSphereRefineAtStart)

  /* ----------------------------------------------------------------------- */
 
  /* set up field names and units -- NOTE: these absolutely MUST be in 
     the same order that they are in Grid_InitializeUniformGrids.C, or 
     else you'll find out that data gets written into incorrectly-named
     fields.  Just FYI. */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;

  if(MetaData.TopGridRank > 1)
    DataLabel[i++] = Vel2Name;

  if(MetaData.TopGridRank > 2)
    DataLabel[i++] = Vel3Name;

  if (TestProblemData.MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = MetalName;

  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RotatingSphereCentralDensity         = %"FSYM"\n"  , RotatingSphereCentralDensity);
    fprintf(Outfptr, "RotatingSphereCentralTemperature         = %"FSYM"\n"  , RotatingSphereCentralTemperature);
    fprintf(Outfptr, "RotatingSphereLambda         = %"FSYM"\n"  , RotatingSphereLambda);
    fprintf(Outfptr, "RotatingSphereTotalEnergy         = %"FSYM"\n"  , RotatingSphereTotalEnergy);
    fprintf(Outfptr, "RotatingSphereCoreRadius         = %"PSYM"\n"  , RotatingSphereCoreRadius);
    fprintf(Outfptr, "RotatingSphereCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RotatingSphereCenterPosition, RotatingSphereCenterPosition+1,
		  RotatingSphereCenterPosition+2);
    fprintf(Outfptr, "RotatingSphereRefineAtStart = %"ISYM"\n", RotatingSphereRefineAtStart);

    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);

    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);

    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);


  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


  if(debug){

    printf("Exiting RotatingSphereInitialize\n");
    fflush(stdout);
  }
 
  return SUCCESS;
 
}
