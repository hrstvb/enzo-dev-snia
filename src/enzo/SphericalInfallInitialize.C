/***********************************************************************
/
/  INITIALIZE A THE SPHERICAL INFALL TEST
/
/  written by: Greg Bryan
/  date:       August, 1995
/  modified1:
/
/  PURPOSE:
/     This routines sets up a test of Bertschinger's 1985 self-similar
/     spherical infall solution (ApJ 1985).
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
#define DEFINE_STORAGE
#include "SphericalInfall.h"
#undef DEFINE_STORAGE
 
/* Set default parameter values. */

const char config_spherical_infall_defaults[] =
"### SPHERICAL INFALL DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    SphericalInfall: {\n"
"        InitialPerturbation 	= 0.1;		# density of central peak\n"
"        SubgridLeft		= 0.0;  	# start of subgrid\n"
"        SubgridRight		= 0.0;		# end of subgrid\n"
"        OmegaBaryonNow 	= 1.0;\n"
"        OmegaCDMNow 		= 0.0;\n"
"        UseBaryons		= TRUE;\n"
"\n"
"        SubgridIsStatic	= FALSE;\n"
"        FixedAcceleration	= FALSE;\n"
"        FixedMass		= -999999.9;\n"
"    };\n"
"};\n";

 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
 
 
int SphericalInfallInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char  line[MAX_LINE_LENGTH];
  int   dim, ret, i;
  int   NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* Error check. */
 
  if (!SelfGravity)
    fprintf(stderr, "SphericalInfall: gravity is off!?!");
 
  /* set default parameters */
 
  float SphericalInfallInitialPerturbation;  	// density of central peak
  FLOAT SphericalInfallSubgridLeft;  		// start of subgrid
  FLOAT SphericalInfallSubgridRight;  		// end of subgrid
  float SphericalInfallOmegaBaryonNow;
  float SphericalInfallOmegaCDMNow;
  int   SphericalInfallUseBaryons;
  int   SphericalInfallSubgridIsStatic;
        SphericalInfallFixedAcceleration;
        SphericalInfallFixedMass;

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    SphericalInfallCenter[dim] = FLOAT_UNDEFINED;

  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_spherical_infall_defaults);


  /* read parameters */
 
  Param.GetScalar(SphericalInfallInitialPerturbation,
		"Problem.SphericalInfall.InitialPerturbation");
  Param.GetScalar(SphericalInfallSubgridLeft,
		"Problem.SphericalInfall.SubgridLeft");
  Param.GetScalar(SphericalInfallOmegaBaryonNow,
		"Problem.SphericalInfall.OmegaBaryonNow");
  Param.GetScalar(SphericalInfallOmegaCDMNow,
		"Problem.SphericalInfall.OmegaCDMNow");
  Param.GetScalar(SphericalInfallSubgridRight,
		"Problem.SphericalInfall.SubgridRight");
  Param.GetScalar(SphericalInfallUseBaryons,
		"Problem.SphericalInfall.UseBaryons");
  Param.GetScalar(SphericalInfallSubgridIsStatic,
		"Problem.SphericalInfall.SubgridIsStatic");
  Param.GetScalar(SphericalInfallFixedAcceleration,
		"Problem.SphericalInfall.FixedAcceleration");
  Param.GetScalar(SphericalInfallFixedMass,
		"Problem.SphericalInfall.FixedMass");
  Param.GetArray(SphericalInfallCenter,
		"Problem.SphericalInfall.Center");
 
 
  /* Set SphericalInfallCenter if left unset. */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    if (SphericalInfallCenter[dim] == FLOAT_UNDEFINED)
      SphericalInfallCenter[dim] =
	(FLOAT(MetaData.TopGridDims[dim]/2) + 0.5) /
	  FLOAT(MetaData.TopGridDims[dim]);
 
  /* If using FixedAcceleration, set FixedMass and zero the Perturbation. */
 
  if (SphericalInfallFixedAcceleration) {
 
    if (SphericalInfallFixedMass == FLOAT_UNDEFINED)
      SphericalInfallFixedMass = SphericalInfallInitialPerturbation /
	POW(float(MetaData.TopGridDims[0]), MetaData.TopGridRank);
 
    SphericalInfallInitialPerturbation = 0.0;
 
  }
 
  /* If there is a subgrid, decrease the perturbation by RefineBy^Rank. */
 
  float TempPerturbation = SphericalInfallInitialPerturbation;
  if (SphericalInfallSubgridLeft != 0.0 || SphericalInfallSubgridRight != 0.0)
    TempPerturbation /= POW(float(RefineBy), MetaData.TopGridRank);
 
  /* Set up grid. */
 
  if (TopGrid.GridData->SphericalInfallInitializeGrid(
					 TempPerturbation,
					 SphericalInfallUseBaryons,
					 SphericalInfallOmegaBaryonNow,
					 SphericalInfallOmegaCDMNow,
					 SphericalInfallSubgridIsStatic
						  ) == FAIL){
    ENZO_FAIL("Error in SphericalInfallInitializeGrid.\n");
  }
 
  /* If requested, create a subgrid */
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((SphericalInfallSubgridRight - SphericalInfallSubgridLeft)/
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
 
    /* compute the dimensions and left/right edges for the subgrid */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
      LeftEdge[dim]    = SphericalInfallSubgridLeft;
      RightEdge[dim]   = SphericalInfallSubgridRight;
    }
 
    /* Set static refine regions. */
 
    if (SphericalInfallSubgridIsStatic == TRUE) {
      StaticRefineRegionLevel[0] = 0;
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	StaticRefineRegionLeftEdge[0][dim] = LeftEdge[dim];
	StaticRefineRegionRightEdge[0][dim] = RightEdge[dim];
      }
      for (dim = MetaData.TopGridRank; dim < MAX_DIMENSION; dim++) {
	StaticRefineRegionLeftEdge[0][dim] = DomainLeftEdge[dim];
	StaticRefineRegionRightEdge[0][dim] = DomainRightEdge[dim];
      }
    }
 
    /* create a new subgrid and initialize it */
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->SphericalInfallInitializeGrid(
                                   SphericalInfallInitialPerturbation,
	                           SphericalInfallUseBaryons,
				   SphericalInfallOmegaBaryonNow,
				   SphericalInfallOmegaCDMNow,
				   FALSE)
	== FAIL) {
      ENZO_FAIL("Error in SphericalInfallInitializeGrid (2).\n");
    }			
  }
 
  /* -------------------------------------------------------------------- */
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set).
     Note: multiply MinimumMassForRefinement by the OmegaBaryonNow since the
     routine that uses this parameter only counts baryonic mass. */
 
  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
 
      MinimumMassForRefinement[i] = SphericalInfallOmegaBaryonNow/
	                            OmegaMatterNow;
      if (CellFlaggingMethod[i] == 4)
	MinimumMassForRefinement[i] = SphericalInfallOmegaCDMNow/
	                              OmegaMatterNow;
 
      MinimumMassForRefinement[i] *= MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *=
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }
  }
 
  /* set up field names and units */
 
  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  DataUnits[5] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "SphericalInfallInitialPerturbation = %"GSYM"\n",
	    SphericalInfallInitialPerturbation);
    fprintf(Outfptr, "SphericalInfallSubgridLeft         = %"GOUTSYM"\n",
	    SphericalInfallSubgridLeft);
    fprintf(Outfptr, "SphericalInfallSubgridRight        = %"GOUTSYM"\n",
	    SphericalInfallSubgridRight);
    fprintf(Outfptr, "SphericalInfallOmegaBaryonNow      = %"GSYM"\n",
	    SphericalInfallOmegaBaryonNow);
    fprintf(Outfptr, "SphericalInfallOmegaCDMNow         = %"GSYM"\n",
	    SphericalInfallOmegaCDMNow);
    fprintf(Outfptr, "SphericalInfallUseBaryons          = %"ISYM"\n",
	    SphericalInfallUseBaryons);
    fprintf(Outfptr, "SphericalInfallSubgridIsStatic     = %"ISYM"\n",
	    SphericalInfallSubgridIsStatic);
    fprintf(Outfptr, "SphericalInfallFixedAcceleration   = %"ISYM"\n",
	    SphericalInfallFixedAcceleration);
    fprintf(Outfptr, "SphericalInfallCenter              = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, SphericalInfallCenter);
    fprintf(Outfptr, "SphericalInfallFixedMass           = %"GSYM"\n\n",
	    SphericalInfallFixedMass);
  }
 
  return SUCCESS;
 
}
