/***********************************************************************
/
/  INITIALIZE A GALAXY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, March 2004
/
/  PURPOSE:
/
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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

const char config_galaxy_simulation_defaults[] =
"### GALAXY SIMULATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    GalaxySimulation: {\n"
"        RefineAtStart 		= TRUE;\n"
"        InitialRefinementLevel	= 0;\n"
"        UseMetallicityField  	= FALSE;\n"
"        InitialTemperature 	= 1000.0;\n"
"        Disk: {\n"
"            Position		= [0.5,0.5,0.5];\n"
"            Radius 		= 0.2;     // [Mpc]\n"
"            Temperature 	= 1.e4;    // [K]\n"
"            ScaleHeightz   	= 325e-6;\n"
"            ScaleHeightR   	= 3500e-6;\n"
"            GasMass		= 4.0e10;\n"
"            GalaxyMass		= 1.0e12;\n"
"        };\n"
"        DarkMatterConcentrationParameter = 12;\n"
"        InflowTime 		= -1;\n"
"        InflowDensity		= 0;\n"
"        AngularMomentum	= [0.0,0.0,0.0];\n"
"        UniformVelocity        = [0.0,0.0,0.0];\n"
"        UniformDensity		= 1.0;\n"
"        UniformEnergy		= 1.0;\n"
"        SubgridLeft		= 0.0;    // start of subgrid(s)\n"
"        SubgridRight		= 0.0;    // end of subgrid(s)\n"
"        
"    };\n"
"};\n";


// Function Prototypes

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);


int GalaxySimulationInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *MetalName = "Metal_Density";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, level, disk, i;

  /* make sure it is 3D */
  
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do GalaxySimulation in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */

  float GalaxySimulationGasMass,
    GalaxySimulationGalaxyMass,    
    GalaxySimulationDiskTemperature,
    GalaxySimulationAngularMomentum[MAX_DIMENSION],
    GalaxySimulationUniformVelocity[MAX_DIMENSION],
    GalaxySimulationUniformDensity,
    GalaxySimulationUniformEnergy;

  FLOAT GalaxySimulationDiskRadius,
    GalaxySimulationDiskPosition[MAX_DIMENSION],
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR;

  float GalaxySimulationInitialTemperature,
    GalaxySimulationDarkMatterConcentrationParameter,
    GalaxySimulationInflowTime,
    GalaxySimulationInflowDensity;

  int   GalaxySimulationRefineAtStart,
    GalaxySimulationInitialRefinementLevel,
    GalaxySimulationUseMetallicityField;
  
  FLOAT GalaxySimulationSubgridLeft, GalaxySimulationSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  int NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

    /* read input from file */

    
    Param.GetScalar(GalaxySimulationRefineAtStart, "Problem.GalaxySimulation.RefineAtStart");
    Param.GetScalar(GalaxySimulationInitialRefinementLevel, "Problem.GalaxySimulation.RefinementLevel");

    Param.GetScalar(GaxySimulationUseMetallicityField, "Problem.GalaxySimulation.UseMetallicityField");
    Param.GetScalar(GalaxySimulationInitialTemperature, "Problem.GalaxySimulation.Temperature");
    Param.GetArray(GalaxySimulationUniformVelocity, "Problem.GalaxySimulation.UniformVelocity");

    // Disk Parameters
    Param.GetScalar(GalaxySimulationDiskRadius, "Problem.GalaxySimulation.Disk.Radius");
    Param.GetScalar(GalaxySimulationGalaxyMass, "Problem.GalaxySimulation.Disk.Mass");
    Param.GetScalar(GalaxySimulationGasMass, "Problem.GalaxySimulation.Disk.GasMass");
    Param.GetArray(GalaxySimulationDiskPosition, "Problem.GalaxySimulation.Disk.Position");
    Param.GetScalar(GalaxySimulationDiskScaleHeightz, "Problem.GalaxySimulation.Disk.ScaleHeightz");
    Param.GetScalar(GalaxySimulationDiskScaleHeightR, "Problem.GalaxySimulation.Disk.ScaleHeightR");
    Param.GetScalar(GalaxySimulationDiskTemperature, "Problem.GalaxySimulation.Disk.Temperature");

    Param.GetScalar(GalaxySimulationDarkMatterConcentrationParameter,
					"Problem.GalaxySimulation.DarkMatterConcentrationParameter");
    Param.GetScalar(GalaxySimulationInflowTime, "Problem.GalaxySimulation.InflowTime");
    Param.GetScalar(GalaxySimulationInflowDensity, "Problem.GalaxySimulation.InflowDensity");
    Param.GetArray(GalaxySimulationAngularMomentum, "Problem.GalaxySimulation.AngularMomentum");

    Param.GetScalar(GalaxySimulationSubgridLeft, "Problem.GalaxySimulation.SubgridLeft");
    Param.GetScalar(GalaxySimulationSubgridRight, "Problem.GalaxySimulation.SubgridRight");
    


  // do a couple of quick error checks
  
  if(GalaxySimulationInitialRefinementLevel > MaximumRefinementLevel){
    ENZO_VFAIL("GalaxySimulationInitialRefinementLevel (%"ISYM") > MaximumRefinementLevel (%"ISYM")!\n",GalaxySimulationInitialRefinementLevel,MaximumRefinementLevel)
  }

  if(GalaxySimulationSubgridLeft > GalaxySimulationSubgridRight){
    ENZO_VFAIL("GalaxySimulationSubgridLeft (%"GOUTSYM") > GalaxySimulationSubgridRight (%"GOUTSYM")!\n")
  }
 
  // Initialize the top grid
  if (TopGrid.GridData->InitializeUniformGrid(GalaxySimulationUniformDensity,
					      GalaxySimulationUniformEnergy,
					      GalaxySimulationUniformEnergy,
					      GalaxySimulationUniformVelocity,
                          ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  /* Create as many subgrids as there are refinement levels
     needed to resolve the initial explosion region upon the start-up. */
 
  HierarchyEntry ** Subgrid;
  if ((MaximumRefinementLevel > 0) && GalaxySimulationRefineAtStart)
    Subgrid   = new HierarchyEntry*[GalaxySimulationInitialRefinementLevel];
  
  /* Create new HierarchyEntries. */
 
  int lev;
  if(GalaxySimulationRefineAtStart)
    for (lev = 0; lev < GalaxySimulationInitialRefinementLevel; lev++)
      Subgrid[lev] = new HierarchyEntry;
  
  if(GalaxySimulationRefineAtStart)
  for (lev = 0; lev < GalaxySimulationInitialRefinementLevel; lev++) {
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((GalaxySimulationSubgridRight - GalaxySimulationSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("GalaxySimulation:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
	     NumberOfSubgridZones[0]);
 
    if (NumberOfSubgridZones[0] > 0) {

      /* fill them out */
 
      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == GalaxySimulationInitialRefinementLevel-1)
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
	LeftEdge[dim]    = GalaxySimulationSubgridLeft;
	RightEdge[dim]   = GalaxySimulationSubgridRight;
      }

      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(GalaxySimulationUniformDensity,
                           GalaxySimulationUniformEnergy,
                           GalaxySimulationUniformEnergy,
                           GalaxySimulationUniformVelocity,
                           ZeroBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }

      /* set up the initial galaxy area on the finest resolution subgrid */
 
      if (lev == GalaxySimulationInitialRefinementLevel - 1)
	if (Subgrid[lev]->GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
						       GalaxySimulationGalaxyMass, 
						       GalaxySimulationGasMass,
						       GalaxySimulationDiskPosition, 
						       GalaxySimulationDiskScaleHeightz,
						       GalaxySimulationDiskScaleHeightR, 
						       GalaxySimulationDarkMatterConcentrationParameter,
						       GalaxySimulationDiskTemperature, 
						       GalaxySimulationInitialTemperature,
						       GalaxySimulationAngularMomentum,
						       GalaxySimulationUniformVelocity,
						       GalaxySimulationUseMetallicityField,
						       GalaxySimulationInflowTime,
						       GalaxySimulationInflowDensity,0)
	    == FAIL) {
	  	  ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
	}
    }
    else{
      printf("GalaxySimulation: single grid start-up.\n");
    }
  }

  /* set up subgrids from level 1 to max refinement level -1 */
 
  for (lev = GalaxySimulationInitialRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  /* set up the root grid */
 
  if (GalaxySimulationInitialRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (TopGrid.GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
							 GalaxySimulationGalaxyMass, 
							 GalaxySimulationGasMass,
							 GalaxySimulationDiskPosition, 
							 GalaxySimulationDiskScaleHeightz,
							 GalaxySimulationDiskScaleHeightR, 
							 GalaxySimulationDarkMatterConcentrationParameter,
							 GalaxySimulationDiskTemperature, 
							 GalaxySimulationInitialTemperature,
							 GalaxySimulationAngularMomentum,
							 GalaxySimulationUniformVelocity,
							 GalaxySimulationUseMetallicityField,
							 GalaxySimulationInflowTime,
							 GalaxySimulationInflowDensity,0) == FAIL) {
            ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
    }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

 /* set up field names and units */

 int count = 0;
 DataLabel[count++] = DensName;
 DataLabel[count++] = TEName;
 if (DualEnergyFormalism)
   DataLabel[count++] = GEName;
 DataLabel[count++] = Vel1Name;
 if(MetaData.TopGridRank > 1)
   DataLabel[count++] = Vel2Name;
 if(MetaData.TopGridRank > 2)
   DataLabel[count++] = Vel3Name;
 if (GalaxySimulationUseMetallicityField)
   DataLabel[count++] = MetalName;

 for (i = 0; i < count; i++)
   DataUnits[i] = NULL;

 /* Write parameters to parameter output file */

 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "GalaxySimulationRefineAtStart      = %"ISYM"\n",
	   GalaxySimulationRefineAtStart);
   fprintf(Outfptr, "GalaxySimulationInitialRefinementLevel      = %"ISYM"\n",
	   GalaxySimulationInitialRefinementLevel);
   fprintf(Outfptr, "GalaxySimulationUseMetallicityField          = %"ISYM"\n",
	   GalaxySimulationUseMetallicityField);
   fprintf(Outfptr, "GalaxySimulationInitialTemperature = %"GOUTSYM"\n",
	   GalaxySimulationInitialTemperature);
   fprintf(Outfptr, "GalaxySimulationUniformVelocity    = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   GalaxySimulationUniformVelocity[0], GalaxySimulationUniformVelocity[1],
	   GalaxySimulationUniformVelocity[2]);
   fprintf(Outfptr, "GalaxySimulationDiskRadius = %"GOUTSYM"\n",
	   GalaxySimulationDiskRadius);
   fprintf(Outfptr, "GalaxySimulationGalaxyMass = %"GOUTSYM"\n",
	   GalaxySimulationGalaxyMass);
   fprintf(Outfptr, "GalaxySimulationGasMass = %"GOUTSYM"\n",
	   GalaxySimulationGasMass);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightz = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightz);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightR = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightR);
   fprintf(Outfptr, "GalaxySimulationDarkMatterConcentrationParameter = %"GOUTSYM"\n",
	   GalaxySimulationDarkMatterConcentrationParameter);
   fprintf(Outfptr, "GalaxySimulationDiskTemperature = %"GOUTSYM"\n",
	   GalaxySimulationDiskTemperature);
   fprintf(Outfptr, "GalaxySimulationInflowTime = %"GOUTSYM"\n",
	   GalaxySimulationInflowTime);
   fprintf(Outfptr, "GalaxySimulationInflowDensity = %"GOUTSYM"\n",
	   GalaxySimulationInflowDensity);
   fprintf(Outfptr, "GalaxySimulationDiskPosition = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationDiskPosition);
   fprintf(Outfptr, "GalaxySimulationAngularMomentum = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationAngularMomentum);
 }

#ifdef USE_MPI

 // BWO: this forces the synchronization of the various point source gravity
 // parameters between processors.  If this is not done, things go to pieces!

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
 MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);

#endif

 return SUCCESS;

}
