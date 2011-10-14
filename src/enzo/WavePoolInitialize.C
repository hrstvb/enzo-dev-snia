/***********************************************************************
/
/  INITIALIZE A WAVE POOL SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/    The wave pool sets up a system which allows a 1d sinusoidal wave to
/    enter from the left boundary.  The initial active region
/    is completely uniform, and wave enters via inflow boundary conditions.
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
#define DEFINE_STORAGE
#include "WavePoolGlobalData.h"
#undef DEFINE_STORAGE

/* Set default parameter values. */

const char config_wave_pool_defaults[] =
"### WAVE POOL DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    WavePool: {\n"
"        Amplitude     = 0.01;   // linear wave\n"
"        Wavelength    = 0.1;    // one-tenth of the box\n"
"        NumberOfWaves = 1;      // just one wave\n"
"        Angle         = 0.0;    // direction of wave propogation wrt x-axis\n"
"\n"
"        Density       = 1.0;    // uniform pool\n"
"        Pressure      = 1.0;"
"        Velocity      = [0.0,0.0,0.0];  // x, y and z velocities
"\n"
"        SubgridLeft   = 0.0;    // start of subgrid\n"
"        SubgridRight  = 0.0;    // end of subgrid\n"
"    };\n"
"};\n";

 
int WavePoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION];
  float WavePoolTotalEnergy;
  FLOAT WavePoolSubgridLeft, WavePoolSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};
 
  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_wave_pool_defaults);

 
  /* read parameters */
 
  Param.GetScalar(WavePoolAmplitude,
		"Problem.WavePool.Amplitude");
  Param.GetScalar(WavePoolWavelength,
                "Problem.WavePool.Wavelength");
  Param.GetScalar(WavePoolNumberOfWaves,
                "Problem.WavePool.NumberOfWaves");
  Param.GetScalar(WavePoolAngle,
                "Problem.WavePool.Angle");
 
  Param.GetScalar(WavePoolDensity,
                "Problem.WavePool.Density");
  Param.GetScalar(WavePoolPressure,
                "Problem.WavePool.Pressure");
  Param.GetArray(WavePoolVelocity,
                "Problem.WavePool.Velocity");
  Param.GetScalar(WavePoolSubgridLeft,
		"Problem.WavePool.SubgridLeft");
  Param.GetScalar(WavePoolSubgridRight,
                "Problem.WavePool.SubgridRight");
 
  /* set the inflow boundary on the left, otherwise leave things alone. */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.LeftFaceBoundaryCondition[dim] = inflow;
 
  /* compute total energy */
 
  WavePoolTotalEnergy = WavePoolPressure/((Gamma - 1.0)*WavePoolDensity);
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    WavePoolTotalEnergy += 0.5*WavePoolVelocity[dim]*WavePoolVelocity[dim];
 
  /* set up grid */
 
  if (TopGrid.GridData->InitializeUniformGrid(WavePoolDensity,
					      WavePoolTotalEnergy,
					      WavePoolTotalEnergy,
					      WavePoolVelocity,
                          ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((WavePoolSubgridRight - WavePoolSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    FLOAT(MetaData.TopGridDims[dim])))
	*RefineBy;
 
  if (NumberOfSubgridZones[0] > 0) {
 
    /* create a new HierarchyEntry, attach to the top grid and fill it out */
 
    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;
 
    /* Compute the dimensions and left/right edges for the subgrid. */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
      LeftEdge[dim]    = WavePoolSubgridLeft;
      RightEdge[dim]   = WavePoolSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->InitializeUniformGrid(WavePoolDensity,
						 WavePoolTotalEnergy,
						 WavePoolTotalEnergy,
						 WavePoolVelocity,
                         ZeroBField) == FAIL) {
            ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
    }			
  }
 
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
    fprintf(Outfptr, "WavePoolAmplitude     = %"FSYM"\n", WavePoolAmplitude);
    fprintf(Outfptr, "WavePoolWavelength    = %"FSYM"\n", WavePoolWavelength);
    fprintf(Outfptr, "WavePoolNumberOfWaves = %"FSYM"\n", WavePoolNumberOfWaves);
    fprintf(Outfptr, "WavePoolAngle         = %"FSYM"\n\n", WavePoolAngle);
 
    fprintf(Outfptr, "WavePoolDensity       = %"FSYM"\n", WavePoolDensity);
    fprintf(Outfptr, "WavePoolPressure      = %"FSYM"\n", WavePoolPressure);
    fprintf(Outfptr, "WavePoolVelocity1     = %"FSYM"\n", WavePoolVelocity[0]);
    fprintf(Outfptr, "WavePoolVelocity2     = %"FSYM"\n", WavePoolVelocity[1]);
    fprintf(Outfptr, "WavePoolVelocity3     = %"FSYM"\n\n", WavePoolVelocity[2]);
 
    fprintf(Outfptr, "WavePoolSubgridLeft   = %"GOUTSYM"\n", WavePoolSubgridLeft);
    fprintf(Outfptr, "WavePoolSubgridRight  = %"GOUTSYM"\n\n", WavePoolSubgridRight);
  }
 
  return SUCCESS;
 
}
