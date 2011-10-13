/***********************************************************************
/
/  INITIALIZE A SHOCK IN A BOX
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
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

const char config_shock_in_a_box_defaults[] = 
"### SHOCK IN A BOX DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    ShockInABox: {\n"
"        Direction = 0;\n"
"        Boundary = 0.5;\n"
"        Density = [-99999.0, -99999.0];\n"
"        Pressure = [-99999.0, -99999.0];\n"
"        Velocity = [-99999.0, -99999.0];\n"
"        SubgridLeft = 0.0;"
"        SubgridRight = 0.0;"
"    };\n"
"};\n";

int ShockInABoxInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION], ShockInABoxDirection;
  float ShockInABoxDensity[2], ShockInABoxPressure[2], ShockInABoxVelocity[2];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* set default parameters */
 
  FLOAT ShockInABoxBoundary = 0.5;     //
 
 
  FLOAT ShockInABoxSubgridLeft  = 0.0;
  FLOAT ShockInABoxSubgridRight = 0.0;
 

  float d1 = 1, v1 = 0, m = 2, p1 = 1, d2, p2, v2, c1, shockspeed = 0;
 
  d2 = d1*((Gamma+1)*m*m)/((Gamma-1)*m*m + 2);
  p2 = p1*(2.0*Gamma*m*m - (Gamma-1))/(Gamma+1);
  c1 = sqrt(Gamma*p1/d1);
  v2 = m*c1*(1-d1/d2);

  shockspeed = 0.9*c1 * m;

  ShockInABoxPressure[0] = p1;
  ShockInABoxPressure[1] = p2;

  ShockInABoxDensity[0] = d1;
  ShockInABoxDensity[1] = d2;

  ShockInABoxVelocity[0] = shockspeed - v1;
  ShockInABoxVelocity[1] = shockspeed - v2;

  Param.SetArray("Problem.ShockInABox.Density",ShockInABoxDensity);
  
  Param.SetArray("Problem.ShockInABox.Pressure", ShockInABoxPressure);
  
  Param.SetArray("Problem.ShockInABox.Velocity", ShockInABoxVelocity);



  /* read parameters */
  
  Param.GetScalar(ShockInABoxDirection, "Problem.ShockInABox.Direction");
  
  Param.GetScalar(ShockInABoxBoundary, "Problem.ShockInABox.Boundary");
  
  Param.GetArray(ShockInABoxDensity, "Problem.ShockInABox.Density");
  
  Param.GetArray(ShockInABoxPressure, "Problem.ShockInABox.Pressure");
  
  Param.GetArray(ShockInABoxVelocity, "Problem.ShockInABox.Velocity");
  
  Param.GetScalar(ShockInABoxSubgridLeft, "Problem.ShockInABox.SubgridLeft");
  Param.GetScalar(ShockInABoxSubgridRight, "Problem.ShockInABox.SubgridRight");


 
 
 
  /* set up grid */
 
  if( ShockInABoxDirection != 0 )
    ENZO_FAIL("Only ShockInABoxDirection=0 supported at the moment!");

//   if (TopGrid.GridData->ShockTubeInitializeGrid(ShockInABoxDirection,
// 						ShockInABoxBoundary,
// 						ShockInABoxDensity,
// 						ShockInABoxPressure,
// 						ShockInABoxVelocity) == FAIL) {
//     ENZO_FAIL("Error in ShockTubeInitializeGrid.\n");
//   }

  if (TopGrid.GridData->
      HydroShockTubesInitializeGrid(ShockInABoxBoundary,
				    ShockInABoxDensity[0], ShockInABoxDensity[1],
				    ShockInABoxVelocity[0], ShockInABoxVelocity[1],
				    0.0, 0.0,
				    0.0, 0.0, 
				    ShockInABoxPressure[0], ShockInABoxPressure[1]) == FAIL) {
    ENZO_FAIL("Error in HydroShockTubesInitializeGrid (called from ShockInABoxInitialize).\n");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((ShockInABoxSubgridRight - ShockInABoxSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    float(MetaData.TopGridDims[dim])))
	*RefineBy;
 
  if (NumberOfSubgridZones[0] > 0) {
 
    /* create a new HierarchyEntry, attach to the top grid and fill it out */
 
    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;
 
    /* compute the dimensions and left/right edges for the subgrid */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
      LeftEdge[dim]    = ShockInABoxSubgridLeft;
      RightEdge[dim]   = ShockInABoxSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->
	HydroShockTubesInitializeGrid(ShockInABoxBoundary,
				      ShockInABoxDensity[0], ShockInABoxDensity[1],
				      ShockInABoxVelocity[0], ShockInABoxVelocity[1],
				      0.0, 0.0,
				      0.0, 0.0, 
				      ShockInABoxPressure[0], ShockInABoxPressure[1]) == FAIL) {
      ENZO_FAIL("Error in HydroShockTubesInitializeGrid (called from ShockInABoxInitialize).\n");
    }
  }
  
  /* Initialize the exterior. */
 
  Exterior.Prepare(TopGrid.GridData);
 
  float InflowValue[5], Dummy[5];
  InflowValue[0] = ShockInABoxDensity[0];
  InflowValue[1] = ShockInABoxPressure[0]/(Gamma-1.0)/ShockInABoxDensity[0]
                   + 0.5*POW(ShockInABoxVelocity[0], 2);
  InflowValue[2] = ShockInABoxVelocity[0];
  InflowValue[3] = 0.0;
  InflowValue[4] = 0.0;
 
  if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
					      Dummy) == FAIL) {
      ENZO_FAIL("Error in InitializeExternalBoundaryFace.\n");
    }
 
  if (MetaData.TopGridRank > 1)
    Exterior.InitializeExternalBoundaryFace(1, reflecting, reflecting,
					    Dummy, Dummy);
  if (MetaData.TopGridRank > 2)
    Exterior.InitializeExternalBoundaryFace(2, reflecting, reflecting,
					    Dummy, Dummy);
 
 
  /* set up field names and units */
 
  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "ShockInABoxDirection     = %"ISYM"\n", ShockInABoxDirection);
    fprintf(Outfptr, "ShockInABoxBoundary      = %"GOUTSYM"\n\n",
	    ShockInABoxBoundary);
 
    fprintf(Outfptr, "ShockInABoxLeftDensity   = %"FSYM"\n", ShockInABoxDensity[0]);
    fprintf(Outfptr, "ShockInABoxLeftPressure  = %"FSYM"\n",
	    ShockInABoxPressure[0]);
    fprintf(Outfptr, "ShockInABoxLeftVelocity  = %"FSYM"\n\n",
	    ShockInABoxVelocity[0]);
 
    fprintf(Outfptr, "ShockInABoxRightDensity  = %"FSYM"\n", ShockInABoxDensity[1]);
    fprintf(Outfptr, "ShockInABoxRightPressure = %"FSYM"\n",
	    ShockInABoxPressure[1]);
    fprintf(Outfptr, "ShockInABoxRightVelocity = %"FSYM"\n\n",
	    ShockInABoxVelocity[1]);
  }
 
  return SUCCESS;
}
