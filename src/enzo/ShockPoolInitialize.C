/***********************************************************************
/
/  INITIALIZE A SHOCK POOL SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/    The shock pool sets up a system which introduces a shock from the
/    left boundary.  The initial active region is completely uniform,
/    and wave enters via inflow boundary conditions.
/
/    In the frame in which the shock is stationary, the definitions are:
/                     |
/       Velocity2     |   Velocity1
/       Density2      |   Density1
/       Pressure2     |   Pressure1
/                     |
/
/    The laboratory frame (in which LabVelocity1 = 0) is related to this
/      frame through:
/                       Velocity1 = LabVelocity1 - ShockVelocity
/                       Velocity2 = LabVelocity2 - ShockVelocity
/
/    The MachNumber = |     Velocity1 / SoundSpeed1 |
/                   = | ShockVelocity / SoundSpeed1 |
/                     (if LabVelocity1 = 0)
/
/    See also Mihalas & Mihalas (Foundations of Radiation Hydrodynamics,
/        p. 236) eq. 56-40
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
#define DEFINE_STORAGE
#include "ShockPoolGlobalData.h"
#undef DEFINE_STORAGE
 
/* Set default parameter values. */

const char config_shock_pool_defaults[] = 
"### SHOCK POOL INITIALIZATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    ShockPool: {\n"
"        Angle         = 0.0;  # x direction\n"
"        MachNumber    = 2.0;  # Velocity[0] / SoundSpeed\n"
"        Density       = 1.0;  # Density in region 1 (preshock)\n"
"        Pressure      = 1.0;  # Pressure in region 1\n"
"        Velocity      = [0.0, 0.0, 0.0];  # LabVelocity in region (these should all be zero for the MachNumber to be correct)\n"
"        SubgridLeft   = 0.0;  # start of subgrid\n"
"        SubgridRight  = 0.0;  # end of subgrid\n"
"    }\n"
"}\n";  

 
int ShockPoolInitialize(FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* parameter declarations */
 
  float ShockPoolMachNumber;
  FLOAT ShockPoolSubgridLeft, ShockPoolSubgridRight;
 
  /* local declarations */
 
  int  dim, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float MachSquared, SoundSpeed1;
  float ShockPoolShockVel, ShockPoolPressure, ShockPoolShockPressure;
  const float TwoPi = 6.283185;
  float ZeroBField[3] = {0.0, 0.0, 0.0} ;

 
  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_shock_pool_defaults);
  
  /* read parameters */
 
  Param.GetScalar(ShockPoolAngle, "Problem.ShockPool.Angle");
  Param.GetScalar(ShockPoolMachNumber, "Problem.ShockPool.MachNumber");
  Param.GetScalar(ShockPoolDensity, "Problem.ShockPool.Density");
  Param.GetScalar(ShockPoolPressure, "Problem.ShockPool.Pressure");

  Param.GetArray(ShockPoolVelocity, "Problem.ShockPool.Velocity");
 
  Param.GetScalar(ShockPoolSubgridLeft, "Problem.ShockPool.SubgridLeft");
  Param.GetScalar(ShockPoolSubgridRight, "Problem.ShockPool.SubgridRight");

 
  /* Compute the physical variables in the postshock region */
 
  MachSquared            = ShockPoolMachNumber * ShockPoolMachNumber;
  ShockPoolShockDensity  = ShockPoolDensity *
                           ((Gamma + 1.0) * MachSquared      ) /
                           ((Gamma - 1.0) * MachSquared + 2.0);
  ShockPoolShockPressure = ShockPoolPressure *
                           (2.0 * Gamma * MachSquared - (Gamma - 1.0)) /
                           (Gamma + 1.0);
  SoundSpeed1 = sqrt(Gamma * ShockPoolPressure / ShockPoolDensity);
  ShockPoolShockVel = SoundSpeed1 * ShockPoolMachNumber *
                      (1.0 - ShockPoolDensity / ShockPoolShockDensity);
  ShockPoolShockVelocity[0] = cos(ShockPoolAngle*TwoPi/360.)*ShockPoolShockVel;
  ShockPoolShockVelocity[1] = sin(ShockPoolAngle*TwoPi/360.)*ShockPoolShockVel;
  ShockPoolShockVelocity[2] = 0.0;
 
  /* Compute total energies */
 
  ShockPoolTotalEnergy = ShockPoolPressure/((Gamma - 1.0)*ShockPoolDensity);
  ShockPoolShockTotalEnergy = ShockPoolShockPressure/((Gamma - 1.0)*
						      ShockPoolShockDensity);
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    ShockPoolTotalEnergy      += 0.5*POW(ShockPoolVelocity[dim], 2);
    ShockPoolShockTotalEnergy += 0.5*POW(ShockPoolShockVelocity[dim], 2);
  }
 
  /* Compute the speed of the shock itself. */
 
  ShockPoolShockSpeed = SoundSpeed1 * ShockPoolMachNumber;
 
  /* set the inflow boundary on the left, otherwise leave things alone. */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.LeftFaceBoundaryCondition[dim] = inflow;
 
  /* set up grid */
 
  if (TopGrid.GridData->InitializeUniformGrid(ShockPoolDensity,
					      ShockPoolTotalEnergy,
					      ShockPoolTotalEnergy,
					      ShockPoolVelocity, ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((ShockPoolSubgridRight - ShockPoolSubgridLeft)/
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
      LeftEdge[dim]    = ShockPoolSubgridLeft;
      RightEdge[dim]   = ShockPoolSubgridRight;
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->InitializeUniformGrid(ShockPoolDensity,
						 ShockPoolTotalEnergy,
						 ShockPoolTotalEnergy,
						 ShockPoolVelocity, ZeroBField) == FAIL) {
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
    fprintf(Outfptr, "ShockPoolAngle        = %"FSYM"\n"  , ShockPoolAngle);
    fprintf(Outfptr, "ShockPoolMachNumber   = %"FSYM"\n\n", ShockPoolMachNumber);
 
    fprintf(Outfptr, "ShockPoolDensity      = %"FSYM"\n"  , ShockPoolDensity);
    fprintf(Outfptr, "ShockPoolPressure     = %"FSYM"\n"  , ShockPoolPressure);
    fprintf(Outfptr, "ShockPoolVelocity1    = %"FSYM"\n"  , ShockPoolVelocity[0]);
    fprintf(Outfptr, "ShockPoolVelocity2    = %"FSYM"\n"  , ShockPoolVelocity[1]);
    fprintf(Outfptr, "ShockPoolVelocity3    = %"FSYM"\n\n", ShockPoolVelocity[2]);
 
    fprintf(Outfptr, "ShockPoolSubgridLeft  = %"GOUTSYM"\n"  , ShockPoolSubgridLeft);
    fprintf(Outfptr, "ShockPoolSubgridRight = %"GOUTSYM"\n\n", ShockPoolSubgridRight);
  }
 
  /* For Zeus solver, subtract kinetic component from TotalEnergy. */
 
  if (HydroMethod == Zeus_Hydro)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      ShockPoolTotalEnergy      -= 0.5*POW(ShockPoolVelocity[dim], 2);
      ShockPoolShockTotalEnergy -= 0.5*POW(ShockPoolShockVelocity[dim], 2);
    }
 
  return SUCCESS;
 
}
