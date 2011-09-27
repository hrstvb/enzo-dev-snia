/***********************************************************************
/
/  INITIALIZE FREE EXPANSION STAGE OF A BLAST WAVE
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:  
/
/  PURPOSE:
/
/   REFERENCE: Draine & Woods, 1991, ApJ
/              Truelove & McKee, 1999, ApJS, 120, 299
/
/   Two dimensional parameters: explosion energy E and ambient density rho_1
/   Two independent variables: radius r, time t
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

/* Set default parameter values. */

const char config_free_expansion_defaults[] = 
"### FREE EXPANSION INITIALIZATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    FreeExpansion: {\n"
"        FullBox      = False;\n"
"        Mass         = 1.0;  # Msun\n"
"        Radius       = 0.1;\n"
"        Density      = 1.0;\n"
"        Energy       = 1e51;  # ergs\n"
"        MaxVelocity  = -99999.0;  # km/s\n"
"        Temperature  = 100;  # K\n"
"        Velocity     = [0.0, 0.0, 0.0];  # initially at rest\n"
"        BField       = [0.0, 0.0, 0.0];  # Gauss\n"
"        SubgridLeft  = 0.0;  # start of subgrid(s)\n"
"        SubgridRight = 0.0;  # end of subgrid(s)\n"
"    }\n"
"}\n";  

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int FreeExpansionInitialize(FILE *Outfptr, HierarchyEntry &TopGrid,
			    TopGridData &MetaData)
{
  char	*DensName  = "Density";
  char	*TEName	   = "TotalEnergy";
  char	*GEName	   = "GasEnergy";
  char	*Vel1Name  = "x-velocity";
  char	*Vel2Name  = "y-velocity";
  char	*Vel3Name  = "z-velocity";
  char	*B1Name	   = "Bx";
  char	*B2Name	   = "By";
  char	*B3Name	   = "Bz";
  char	*PhiName   = "Phi";
  char	*DebugName = "Debug";
  char	*Phi_pName = "Phip";

  /* parameter declarations */

  float FreeExpansionTotalEnergy, FreeExpansionGasEnergy, B2, V2;
  FLOAT FreeExpansionSubgridLeft, FreeExpansionSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  int  dim, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

  int FreeExpansionFullBox;
  float FreeExpansionVelocity[3];
  float FreeExpansionBField[3];
  float FreeExpansionDensity;
  float FreeExpansionMaxVelocity;
  float FreeExpansionRadius;
  float FreeExpansionMass;
  double FreeExpansionEnergy;
  float FreeExpansionTemperature;

  float dx = (DomainRightEdge[0] - DomainLeftEdge[0]) / MetaData.TopGridDims[0];

  /* Use 3.5 zones on the finest level to resolve the initial explosion at t=0. */
  float dr = 3.5*dx*max(POW(RefineBy,-MaximumRefinementLevel), 0.25);
  
  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_free_expansion_defaults);

  /* read parameters */

  Param.GetScalar(FreeExpansionFullBox, "Problem.FreeExpansion.FullBox");
  Param.GetScalar(FreeExpansionMass, "Problem.FreeExpansion.Mass");
  Param.GetScalar(FreeExpansionRadius, "Problem.FreeExpansion.Radius");
  Param.GetScalar(FreeExpansionDensity, "Problem.FreeExpansion.Density");
  Param.GetScalar(FreeExpansionEnergy, "Problem.FreeExpansion.Energy");
  Param.GetScalar(FreeExpansionMaxVelocity, "Problem.FreeExpansion.MaxVelocity");
  Param.GetScalar(FreeExpansionTemperature, "Problem.FreeExpansion.Temperature");
  Param.GetScalar(FreeExpansionSubgridLeft, "Problem.FreeExpansion.SubgridLeft");
  Param.GetScalar(FreeExpansionSubgridRight, "Problem.FreeExpansion.SubgridRight");

  Param.GetArray(FreeExpansionVelocity, "Problem.FreeExpansion.Velocity");
  Param.GetArray(FreeExpansionBField, "Problem.FreeExpansion.BField");


  /* get units */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, PressureUnits, MagneticUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, 0.0);

  PressureUnits = DensityUnits * (LengthUnits/TimeUnits)*(LengthUnits/TimeUnits);
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    FreeExpansionBField[dim] /= MagneticUnits;

  if (debug) 
    printf("Bunits = %"GSYM" G, Bfield(code) = %"GSYM" %"GSYM" %"GSYM"\n",
	   MagneticUnits, FreeExpansionBField[0], FreeExpansionBField[1],
	   FreeExpansionBField[2]);

  /* Set up current problem time, ambient total energy. */

  MetaData.Time         = 0.0;
  FreeExpansionGasEnergy = FreeExpansionTemperature / TemperatureUnits / 
    ((Gamma-1.0)*DEFAULT_MU);

  V2 = 0;
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    V2 += FreeExpansionVelocity[dim] * FreeExpansionVelocity[dim];
  FreeExpansionTotalEnergy = FreeExpansionGasEnergy + 0.5 * V2;

  if (HydroMethod == MHD_RK) {
    B2 = 0.0;
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      B2 += FreeExpansionBField[dim] * FreeExpansionBField[dim];
    FreeExpansionTotalEnergy + 0.5 * B2 / FreeExpansionDensity;
  }

  
  /* set periodic boundaries (if FullBox=1),
     otherwise keep reflecting (the default) */

  if (FreeExpansionFullBox)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
      MetaData.RightFaceBoundaryCondition[dim] = periodic;
    }

  /* set up uniform grid as of before explosion */

  TopGrid.GridData->InitializeUniformGrid(FreeExpansionDensity, 
					  FreeExpansionTotalEnergy,
					  FreeExpansionGasEnergy,
					  FreeExpansionVelocity,
					  FreeExpansionBField);

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
	nint((FreeExpansionSubgridRight - FreeExpansionSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *POW(RefineBy, lev + 1);

    if (debug)
      printf("FreeExpansion:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1, 
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
	LeftEdge[dim]    = FreeExpansionSubgridLeft;
	RightEdge[dim]   = FreeExpansionSubgridRight;
      }

      /* create a new subgrid and initialize it */

      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
					  LeftEdge, RightEdge, 0);
      Subgrid[lev]->GridData->InitializeUniformGrid(FreeExpansionDensity,
						    FreeExpansionTotalEnergy,
						    FreeExpansionGasEnergy,
						    FreeExpansionVelocity,
						    FreeExpansionBField);

      /* set up the initial explosion area on the finest resolution subgrid */

      if (lev == MaximumRefinementLevel - 1)
	Subgrid[lev]->GridData->
	  FreeExpansionInitializeGrid(FreeExpansionFullBox, FreeExpansionDensity,
				      FreeExpansionEnergy, FreeExpansionMaxVelocity,
				      FreeExpansionMass, FreeExpansionRadius, 
				      DensityUnits, VelocityUnits, LengthUnits,
				      TimeUnits);

    } else
      printf("FreeExpansion: single grid start-up.\n");
  } // ENDFOR levels

  /* set up subgrids from level 1 to max refinement level -1 */

  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    
    Subgrid[lev]->GridData->ProjectSolutionToParentGrid(*(Subgrid[lev-1]->GridData));
  
  /* set up the root grid */

  if (MaximumRefinementLevel > 0)
    Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData));
  else
    TopGrid.GridData->
      FreeExpansionInitializeGrid(FreeExpansionFullBox, FreeExpansionDensity,
				  FreeExpansionEnergy, FreeExpansionMaxVelocity,
				  FreeExpansionMass, FreeExpansionRadius, 
				  DensityUnits, VelocityUnits, LengthUnits, 
				  TimeUnits);

  /* set up field names and units */
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = B1Name;
    DataLabel[i++] = B2Name;
    DataLabel[i++] = B3Name;
    DataLabel[i++] = PhiName;
    if (UseDivergenceCleaning) {
      DataLabel[i++] = Phi_pName;
      DataLabel[i++] = DebugName;
    }
  }

  for (int j=0; j< i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "FreeExpansionFullBox         = %"ISYM"\n", FreeExpansionFullBox);
    fprintf(Outfptr, "FreeExpansionDensity         = %"FSYM"\n", FreeExpansionDensity);
    fprintf(Outfptr, "FreeExpansionMass            = %"FSYM"\n", FreeExpansionMass);
    fprintf(Outfptr, "FreeExpansionRadius          = %"FSYM"\n", FreeExpansionRadius);
    fprintf(Outfptr, "FreeExpansionTemperature     = %"GSYM"\n", FreeExpansionTemperature);
    fprintf(Outfptr, "FreeExpansionEnergy          = %lg\n"    , FreeExpansionEnergy);
    fprintf(Outfptr, "FreeExpansionMaxVelocity     = %"FSYM"\n", FreeExpansionMaxVelocity);
  }

  return SUCCESS;

}
