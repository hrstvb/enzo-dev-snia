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
 
// This routine intializes a new simulation based on the parameter file.
//
 
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
 
int RotatingSphereInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
  if(debug){
    printf("Entering RotatingSphereInitialize\n");
    fflush(stdout);
  }

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

  /* parameter declarations */
 
  FLOAT RotatingSphereSubgridLeft[MAX_DIMENSION], RotatingSphereSubgridRight[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT RotatingSphereCenterPosition[MAX_DIMENSION];

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 3D */
 
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do RotatingSphere in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }
 
  for(i=0; i<MAX_DIMENSION; i++)
    RotatingSphereCenterPosition[i] = 0.5;  // right in the middle of the box

  float RotatingSphereVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float RotatingSphereBField[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  FLOAT RotatingSphereRadius = 0.3;
  float RotatingSphereLambda = 0.05;
  float RotatingSphereOverdensity = 20.0;
  float RotatingSphereDensity = 1.0;
  float RotatingSphereTotalEnergy = 1.0;
  float Pi                      = 3.14159;

  /* set no subgrids by default. */
 
  RotatingSphereSubgridLeft[0] = RotatingSphereSubgridLeft[1] = 
    RotatingSphereSubgridLeft[2] = 0.0;    // start of subgrid(s)

  RotatingSphereSubgridRight[0] = RotatingSphereSubgridRight[1] = 
    RotatingSphereSubgridRight[2] = 0.0;    // end of subgrid(s)

  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters specifically for radiating shock problem*/

    ret += sscanf(line, "RotatingSphereOverdensity  = %"FSYM, &RotatingSphereOverdensity);
    ret += sscanf(line, "RotatingSphereSubgridLeft = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereSubgridLeft,RotatingSphereSubgridLeft+1,RotatingSphereSubgridLeft+2);
    ret += sscanf(line, "RotatingSphereSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereSubgridRight,RotatingSphereSubgridRight+1,RotatingSphereSubgridRight+2);
    ret += sscanf(line, "RotatingSphereLambda = %"FSYM,
		        &RotatingSphereLambda);

    ret += sscanf(line, "RotatingSphereTotalEnergy = %"FSYM,
		        &RotatingSphereTotalEnergy);

    ret += sscanf(line, "RotatingSphereRadius = %"PSYM,
		        &RotatingSphereRadius);
    ret += sscanf(line, "RotatingSphereCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		  RotatingSphereCenterPosition, RotatingSphereCenterPosition+1,
		  RotatingSphereCenterPosition+2);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && (strstr(line, "RotatingSphere") || strstr(line, "TestProblem")) &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	 "*** warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file
 
 
  if (TopGrid.GridData->InitializeUniformGrid(RotatingSphereDensity,
					      RotatingSphereTotalEnergy,
					      RotatingSphereTotalEnergy,
					      RotatingSphereVelocity,
					      RotatingSphereBField) == FAIL) {
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
	nint((RotatingSphereSubgridRight[dim] - RotatingSphereSubgridLeft[dim])/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("RotatingSphere:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
	LeftEdge[dim]    = RotatingSphereSubgridLeft[dim];
	RightEdge[dim]   = RotatingSphereSubgridRight[dim];
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(RotatingSphereDensity,
						   RotatingSphereTotalEnergy,
						   RotatingSphereTotalEnergy,
							RotatingSphereVelocity,
							RotatingSphereBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->RotatingSphereInitializeGrid(RotatingSphereRadius,
								   RotatingSphereCenterPosition,
								   RotatingSphereLambda,
								   RotatingSphereOverdensity) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in RotatingSphereInitialize[Sub]Grid.");
	}

    }
    else{
      printf("RotatingSphere: single grid start-up.\n");
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
    if (TopGrid.GridData->RotatingSphereInitializeGrid(RotatingSphereRadius,
							 RotatingSphereCenterPosition,
							 RotatingSphereLambda,
							 RotatingSphereOverdensity) == FAIL) {
            ENZO_FAIL("Error in RotatingSphereInitializeGrid.");
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
    fprintf(Outfptr, "RotatingSphereOverdensity         = %"FSYM"\n"  , RotatingSphereOverdensity);
    fprintf(Outfptr, "RotatingSphereLambda         = %"FSYM"\n"  , RotatingSphereLambda);
    fprintf(Outfptr, "RotatingSphereTotalEnergy         = %"FSYM"\n"  , RotatingSphereTotalEnergy);
    fprintf(Outfptr, "RotatingSphereRadius         = %"PSYM"\n"  , RotatingSphereRadius);
    fprintf(Outfptr, "RotatingSphereCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RotatingSphereCenterPosition, RotatingSphereCenterPosition+1,
		  RotatingSphereCenterPosition+2);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


  if(debug){

    printf("Exiting RotatingSphereInitialize\n");
    fflush(stdout);
  }
 
  return SUCCESS;
 
}
