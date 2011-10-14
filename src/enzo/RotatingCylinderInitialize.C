/***********************************************************************
/
/  Rotating Cylinder Test Problem
/
/  written by: Brian O'Shea
/  date:       February 2008
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
 
int RotatingCylinderInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
  if(debug){
    printf("Entering RotatingCylinderInitialize\n");
    fflush(stdout);
  }

  char *DensName = "Density";
  char *TEName   = "Total_Energy";
  char *GEName   = "Gas_Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *MetalName = "Metal_Density";

  /* parameter declarations */
 

  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];


  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 3D */
 
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do RotatingCylinder in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set no subgrids by default. */
  const char config_rotating_cylinder_defaults[] = 
  "### ROTATING CYLINDER DEFAULTS ###\n"
  "\n"
  "Problem: {\n"
  "    RotatingCylinder: {\n"
  "        Velocity = [0.0, 0.0, 0.0];\n"   // gas initally at rest
  "        BField   = [0.0, 0.0, 0.0];\n"   // gas initally at rest
  "        RotatingCylinderCenterPosition = [0.5, 0.5, 0.5];\n" // right in the middle of the box
  "        RotatingCylinderSubgridLeft=[0.0, 0.0, 0.0];\n"
  "        RotatingCylinderSubgridRight=[0.0, 0.0, 0.0];\n"
  "        Radius = 0.3;\n"
  "        Lambda = 0.05;\n"
  "        Overdensity = 20.0;\n"
  "        Density = 1.0;\n"
  "        TotalEnergy = 1.0;\n"
  "    };\n"
  "};\n";

  FLOAT RotatingCylinderSubgridLeft[MAX_DIMENSION];
  FLOAT RotatingCylinderSubgridRight[MAX_DIMENSION];
  FLOAT RotatingCylinderCenterPosition[MAX_DIMENSION];
  float RotatingCylinderVelocity[MAX_DIMENSION];   // gas initally at rest
  float RotatingCylinderBField[MAX_DIMENSION];   // gas initally at rest
  FLOAT RotatingCylinderRadius;
  float RotatingCylinderLambda;
  float RotatingCylinderOverdensity;
  float RotatingCylinderDensity;
  float RotatingCylinderTotalEnergy;

  Param.Update(config_rotating_cylinder_defaults);

  /* read input from file */
 
  /* read parameters specifically for radiating shock problem*/

  Param.GetArray(RotatingCylinderSubgridLeft, "Problem.RotatingCylinder.SubgridLeft");
  Param.GetArray(RotatingCylinderSubgridRight, "Problem.RotatingCylinder.SubgridRight");
  Param.GetArray(RotatingCylinderCenterPosition, "Problem.RotatingCylinder.CenterPosition");
  Param.GetScalar(RotatingCylinderOverdensity, "Problem.RotatingCylinder.Overdensity");
  Param.GetScalar(RotatingCylinderLambda, "Problem.RotatingCylinder.Lambda");
  Param.GetScalar(RotatingCylinderTotalEnergy, "Problem.RotatingCylinder.TotalEnergy");
  Param.GetScalar(RotatingCylinderRadius, "Problem.RotatingCylinder.Radius");
  Param.GetScalar(TestProblemData.UseMetallicityField, "Problem.RotatingCylinder.UseMetallicityField");
  Param.GetScalar(TestProblemData.MetallicityField_Fraction, "Problem.RotatingCylinder.InitialMetallicityFraction"); 
  if (TopGrid.GridData->InitializeUniformGrid(RotatingCylinderDensity,
					      RotatingCylinderTotalEnergy,
					      RotatingCylinderTotalEnergy,
					      RotatingCylinderVelocity,
					      RotatingCylinderBField) == FAIL){
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* Create as many subgrids as there are refinement levels
     needed to resolve the initial explosion region upon the start-up. */
 
  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0)
    Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
 
  /* Create new HierarchyEntries. */
 
  int lev;
  for (lev = 0; lev < MaximumRefinementLevel; lev++)
    Subgrid[lev] = new HierarchyEntry;
 
  for (lev = 0; lev < MaximumRefinementLevel; lev++) {
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((RotatingCylinderSubgridRight[dim] - RotatingCylinderSubgridLeft[dim])/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("RotatingCylinder:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
	     NumberOfSubgridZones[0]);
 
    if (NumberOfSubgridZones[0] > 0) {
 
      /* fill them out */
 
      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == MaximumRefinementLevel-1)
	Subgrid[lev]->NextGridNextLevel = NULL;
      else
	Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
      if (lev == 0)
	Subgrid[lev]->ParentGrid        = &TopGrid;
      else
	Subgrid[lev]->ParentGrid        = Subgrid[lev-1];
 
      /* compute the dimensions and left/right edges for the subgrid */
 
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
	LeftEdge[dim]    = RotatingCylinderSubgridLeft[dim];
	RightEdge[dim]   = RotatingCylinderSubgridRight[dim];
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(RotatingCylinderDensity,
						   RotatingCylinderTotalEnergy,
						   RotatingCylinderTotalEnergy,
							RotatingCylinderVelocity,
							RotatingCylinderBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->RotatingCylinderInitializeGrid(RotatingCylinderRadius,
								   RotatingCylinderCenterPosition,
								   RotatingCylinderLambda,
								   RotatingCylinderOverdensity) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in RotatingCylinderInitialize[Sub]Grid.");
	}

    }
    else{
      printf("RotatingCylinder: single grid start-up.\n");
    }
  }

 
  /* set up subgrids from level 1 to max refinement level -1 */
 
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
 
  /* set up the root grid */
 
  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (TopGrid.GridData->RotatingCylinderInitializeGrid(RotatingCylinderRadius,
							 RotatingCylinderCenterPosition,
							 RotatingCylinderLambda,
							 RotatingCylinderOverdensity) == FAIL) {
            ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
    }

 
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

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = MetalName;

  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RotatingCylinderOverdensity         = %"FSYM"\n"  , RotatingCylinderOverdensity);
    fprintf(Outfptr, "RotatingCylinderLambda         = %"FSYM"\n"  , RotatingCylinderLambda);
    fprintf(Outfptr, "RotatingCylinderTotalEnergy         = %"FSYM"\n"  , RotatingCylinderTotalEnergy);
    fprintf(Outfptr, "RotatingCylinderRadius         = %"PSYM"\n"  , RotatingCylinderRadius);
    fprintf(Outfptr, "RotatingCylinderCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RotatingCylinderCenterPosition, RotatingCylinderCenterPosition+1,
		  RotatingCylinderCenterPosition+2);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


  if(debug){

    printf("Exiting RotatingCylinderInitialize\n");
    fflush(stdout);
  }
 
  return SUCCESS;
 
}
