/***********************************************************************
/
/  INITIALIZE A COSMOLOGY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/  date:       July 2003
/
/  PURPOSE:  Initialize for cosmology simulations.  Reads in a number
/      of initial grids.  If more than one, then they are all numbered.
/      We currently assume that all are subgrids except the first.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file
 
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
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "fortran.def"
#include "CommunicationUtilities.h"

/* Set default parameter values. */

const char config_nested_cosmology_simulation_defaults[] = 
"### NESTED COSMOLOGY SIMULATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    CosmologySimulation: {\n"
"        OmegaBaryonNow       = 1.0;\n"
"        OmegaCDMNow          = 0.0;\n"
"        InitialTemperature   = -99999.0;\n"
"\n" 
"        DensityName           = \"\";\n"
"        VelocitiesNames       = \"\";\n"
"        TotalEnergyName       = \"\";\n"
"        GasEnergyName         = \"\";\n"
"        ParticlePositionName  = \"\";\n"
"        ParticleVelocityName  = \"\";\n"
"        ParticleDisplacementName = \"\";\n"
"        ParticleMassName      = \"\";\n"
"        ParticleTypeName      = \"\";\n"
"        VelocityNames         = [\"\", \"\", \"\"];\n"
"        ParticleVelocityNames = [\"\", \"\", \"\"];\n"
"        ParticlePositionNames = [\"\", \"\", \"\"];\n"
"        ParticleDisplacementNames = [\"\", \"\", \"\"];\n"
"\n" 
"        InitialFractionHII   = 1.2e-5;\n"
"        InitialFractionHeII  = 1.0e-14;\n"
"        InitialFractionHeIII = 1.0e-17;\n"
"        InitialFractionHM    = 2.0e-9;\n"
"        InitialFractionH2I   = 2.0e-20;\n"
"        InitialFractionH2II  = 3.0e-14;\n"
"        InitialFractionMetal = 1.0e-10;\n"
"        UseMetallicityField  = False;\n"
"\n"  
"        ManuallySetParticleMassRatio = False;\n"
"        ManualParticleMassRatio = 1.0;\n"
"\n"  		
"        CalculatePositions   = False; \n"
"\n"  
"        InitialUniformBField = [0.0, 0.0, 0.0];  # in proper Gauss\n"
"\n"  
"        RadHydroInitialRadiationEnergy = 1.0e-32;\n"
"\n"  
"        SubgridsAreStatic    = True;\n"
"        InitialSubgrids = [];\n"
"    };\n"
"};\n";


// Function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationBroadcastValue(PINT *Value, int BroadcastProcessor);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
// Cosmology Parameters (that need to be shared)
 
static float CosmologySimulationOmegaBaryonNow;
static float CosmologySimulationOmegaCDMNow;
static float CosmologySimulationInitialTemperature;
 
static char *CosmologySimulationDensityName;
static char *CosmologySimulationTotalEnergyName;
static char *CosmologySimulationGasEnergyName;
static char *CosmologySimulationParticlePositionName;
static char *CosmologySimulationParticleVelocityName;
static char *CosmologySimulationParticleDisplacementName;
static char *CosmologySimulationParticleMassName;
static char *CosmologySimulationParticleTypeName;
static char *CosmologySimulationVelocityNames[MAX_DIMENSION];
static char *CosmologySimulationParticleVelocityNames[MAX_DIMENSION];
static char *CosmologySimulationParticlePositionNames[MAX_DIMENSION];
static char *CosmologySimulationParticleDisplacementNames[MAX_DIMENSION];
 
static int   CosmologySimulationSubgridsAreStatic;
 
static float CosmologySimulationInitialFractionHII;
static float CosmologySimulationInitialFractionHeII;
static float CosmologySimulationInitialFractionHeIII;
static float CosmologySimulationInitialFractionHM;
static float CosmologySimulationInitialFractionH2I;
static float CosmologySimulationInitialFractionH2II;
static float CosmologySimulationInitialFractionMetal;
static int   CosmologySimulationUseMetallicityField;
 
static int CosmologySimulationManuallySetParticleMassRatio;
static float CosmologySimulationManualParticleMassRatio;

static int   CosmologySimulationCalculatePositions;

static float CosmologySimulationInitialUniformBField[MAX_DIMENSION];  // in proper Gauss

#define MAX_INITIAL_GRIDS 10
 
static  int   CosmologySimulationGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
static  int   CosmologySimulationGridLevel[MAX_INITIAL_GRIDS];
static  FLOAT CosmologySimulationGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
static  FLOAT CosmologySimulationGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
 

 
int NestedCosmologySimulationInitialize(FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData)
{
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
  char *GPotName  = "Grav_Potential";
  char *MetalName = "Metal_Density";
  char *ForbidName = "ForbiddenRefinement";
  char *MachName   = "Mach";
  char *CRName     = "CR_Density";
  char *PSTempName = "PreShock_Temperature";
  char *PSDenName  = "PreShock_Density";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";

 
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};
 
  // Declarations
 
  char *dummy = new char[MAX_LINE_LENGTH];
  char dummy_arr[MAX_DIMENSION][MAX_LINE_LENGTH];

  int i, j, dim, gridnum, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
 
  char *DensityName = NULL,
    *TotalEnergyName = NULL,
    *GasEnergyName = NULL,
    *ParticlePositionName = NULL,
    *ParticleVelocityName = NULL, 
    *ParticleDisplacementName = NULL,
    *ParticleMassName = NULL, 
    *VelocityNames[MAX_DIMENSION],
    *ParticleTypeName = NULL, 
    *ParticlePositionNames[MAX_DIMENSION],
    *ParticleVelocityNames[MAX_DIMENSION],
    *ParticleDisplacementNames[MAX_DIMENSION];
 
  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_nested_cosmology_simulation_defaults);

  // Set all char arrays to NULL
  CosmologySimulationDensityName          = NULL;
  CosmologySimulationTotalEnergyName      = NULL;
  CosmologySimulationGasEnergyName        = NULL;
  CosmologySimulationParticlePositionName = NULL;
  CosmologySimulationParticleVelocityName = NULL;
  CosmologySimulationParticleDisplacementName = NULL;
  CosmologySimulationParticleMassName     = NULL;
  CosmologySimulationParticleTypeName     = NULL;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    ParticleVelocityNames[dim] = NULL;
    ParticlePositionNames[dim] = NULL;
    ParticleDisplacementNames[dim] = NULL;
    CosmologySimulationParticleVelocityNames[dim] = NULL;
    CosmologySimulationInitialUniformBField[dim] = 0.0;
    CosmologySimulationParticlePositionNames[dim] = NULL;
    CosmologySimulationParticleDisplacementNames[dim] = NULL;
    CosmologySimulationVelocityNames[dim] = NULL;
  }
 
  // Set default parameters: parameters, names and subgrid info 
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    CosmologySimulationGridLevel[i] = 1;
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    CosmologySimulationGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    CosmologySimulationGridRightEdge[0][dim] = DomainRightEdge[dim];
    CosmologySimulationGridDimension[0][dim] = MetaData.TopGridDims[dim];
  } 
  CosmologySimulationGridLevel[0] = 0;

  dummy[0] = 0;
  for( i=0; i<MAX_DIMENSION; i++) dummy_arr[i][0] = 0;
 

  // Error check
 
  if (!ComovingCoordinates) {
    ENZO_FAIL("ComovingCoordinates must be TRUE!\n");
  }
 
  if (DualEnergyFormalism == FALSE && HydroMethod != Zeus_Hydro)
    fprintf(stderr, "CosmologySimulation: DualEnergyFormalism is off!\n");
  if (!SelfGravity)
    fprintf(stderr, "CosmologySimulation: gravity is off!?!\n");
 

  // Read parameters
 
  Param.GetScalar(CosmologySimulationOmegaBaryonNow, "Problem.CosmologySimulation.OmegaBaryonNow");
  Param.GetScalar(CosmologySimulationOmegaCDMNow, "Problem.CosmologySimulation.OmegaCDMNow");
  Param.GetScalar(CosmologySimulationInitialTemperature, "Problem.CosmologySimulation.InitialTemperature");

  Param.GetScalar(dummy, "Problem.CosmologySimulation.DensityName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationDensityName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationDensityName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.VelocitiesName");
  if( strlen(dummy) > 0 ) {
    for( i=0; i<MAX_DIMENSION; i++) {
      CosmologySimulationVelocityNames[i] = new char[MAX_LINE_LENGTH];
      strcpy(CosmologySimulationVelocityNames[i], dummy);
    }
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.TotalEnergyName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationTotalEnergyName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationTotalEnergyName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.GasEnergyName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationGasEnergyName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationGasEnergyName, dummy);
  }
  
  Param.GetScalar(dummy, "Problem.CosmologySimulation.ParticlePositionName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationParticlePositionName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationParticlePositionName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.ParticleVelocityName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationParticleVelocityName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationParticleVelocityName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.ParticleMassName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationParticleMassName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationParticleMassName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.ParticleTypeName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationParticleTypeName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationParticleTypeName, dummy);
  }

  Param.GetScalar(dummy, "Problem.CosmologySimulation.ParticleDisplacementName");
  if( strlen(dummy) > 0 ) {
    CosmologySimulationParticleDisplacementName = new char[MAX_LINE_LENGTH];
    strcpy(CosmologySimulationParticleDisplacementName, dummy);
  }

  Param.GetArray(dummy_arr, "Problem.CosmologySimulation.VelocityNames");  
  for( i=0; i<MAX_DIMENSION; i++) {
    if( strlen(dummy_arr[i]) > 0 ) {
      if( CosmologySimulationVelocityNames[i] == NULL )
	CosmologySimulationVelocityNames[i] = new char[MAX_LINE_LENGTH];
      strcpy(CosmologySimulationVelocityNames[i], dummy_arr[i]);
    }
  }

  Param.GetArray(dummy_arr, "Problem.CosmologySimulation.ParticlePositionNames");
  for( i=0; i<MAX_DIMENSION; i++) {
    if( strlen(dummy_arr[i]) > 0 ) {
      CosmologySimulationParticlePositionNames[i] = new char[MAX_LINE_LENGTH];
      strcpy(CosmologySimulationParticlePositionNames[i], dummy_arr[i]);
    }
  }

  Param.GetArray(dummy_arr, "Problem.CosmologySimulation.VelocityNames");
  for( i=0; i<MAX_DIMENSION; i++) {
    if( strlen(dummy_arr[i]) > 0 ) {
      CosmologySimulationParticleVelocityNames[i] = new char[MAX_LINE_LENGTH];
      strcpy(CosmologySimulationParticleVelocityNames[i], dummy_arr[i]);
    }
  }

  Param.GetArray(dummy_arr, "Problem.CosmologySimulation.ParticleDisplacementNames");
  for( i=0; i<MAX_DIMENSION; i++) {
    if( strlen(dummy_arr[i]) > 0 ) {
      CosmologySimulationParticleDisplacementNames[i] = new char[MAX_LINE_LENGTH];
      strcpy(CosmologySimulationParticleDisplacementNames[i], dummy_arr[i]);
    }
  }

  Param.GetScalar(CosmologySimulationSubgridsAreStatic, "Problem.CosmologySimulation.SubgridsAreStatic");
  
  int NumberOfInitialSubgrids = Param.Size("Problem.CosmologySimulation.InitialSubgrids");
  if (NumberOfInitialSubgrids+1 > MAX_INITIAL_GRIDS) {
    ENZO_VFAIL("You've exceeded the maximum number of CosmologySimulation initial grids (%d)!\n",MAX_INITIAL_GRIDS)
      }
  CosmologySimulationNumberOfInitialGrids = NumberOfInitialSubgrids + 1;  // add 1 for the root grid
  
  char InitialSubgridNames[MAX_LINE_LENGTH][MAX_INITIAL_GRIDS];
  Param.GetArray(InitialSubgridNames,"Problem.CosmologySimulation.InitialSubgrids");
  
  for (i = 0; i < NumberOfInitialSubgrids; i++) {
    Param.GetArray(CosmologySimulationGridLeftEdge[i+1], "Problem.CosmologySimulation.%s.LeftEdge",InitialSubgridNames[i]);
    Param.GetArray(CosmologySimulationGridRightEdge[i+1], "Problem.CosmologySimulation.%s.RightEdge",InitialSubgridNames[i]);
    Param.GetArray(CosmologySimulationGridDimension[i+1], "Problem.CosmologySimulation.%s.Dimension",InitialSubgridNames[i]);
    Param.GetScalar(CosmologySimulationGridLevel[i+1], "Problem.CosmologySimulation.%s.Level",InitialSubgridNames[i]);
  }
  
  
  Param.GetScalar(CosmologySimulationInitialFractionHII, "Problem.CosmologySimulation.InitialFractionHII");
  Param.GetScalar(CosmologySimulationInitialFractionHeII, "Problem.CosmologySimulation.InitialFractionHeII");
  Param.GetScalar(CosmologySimulationInitialFractionHeIII, "Problem.CosmologySimulation.InitialFractionHeIII");
  Param.GetScalar(CosmologySimulationInitialFractionHM, "Problem.CosmologySimulation.InitialFractionHM");
  Param.GetScalar(CosmologySimulationInitialFractionH2I, "Problem.CosmologySimulation.InitialFractionH2I");
  Param.GetScalar(CosmologySimulationInitialFractionH2II, "Problem.CosmologySimulation.InitialFractionH2II");
  Param.GetScalar(CosmologySimulationInitialFractionMetal, "Problem.CosmologySimulation.InitialFractionMetal");
  Param.GetScalar(CosmologySimulationUseMetallicityField, "Problem.CosmologySimulation.UseMetallicityField");
  
  Param.GetScalar(CosmologySimulationManuallySetParticleMassRatio, "Problem.CosmologySimulation.ManuallySetParticleMassRatio");
  Param.GetScalar(CosmologySimulationManualParticleMassRatio, "Problem.CosmologySimulation.ManualParticleMassRatio");
  
  Param.GetScalar(CosmologySimulationCalculatePositions, "Problem.CosmologySimulation.CalculatePositions");
  
  Param.GetArray(CosmologySimulationInitialUniformBField, "Problem.CosmologySimulation.InitialUniformBField");

 
  // More error checking
 
  if (CosmologySimulationDensityName == NULL &&
      (CosmologySimulationParticlePositionName == NULL &&
       CosmologySimulationParticleDisplacementName == NULL &&
       !CosmologySimulationCalculatePositions)) {
    ENZO_FAIL("Missing initial data.\n");
  }
 
  if (CosmologySimulationDensityName != NULL && CellFlaggingMethod[0] != 2)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");
 
  if (CosmologySimulationDensityName == NULL && CellFlaggingMethod[0] != 4)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");
 
  if (CosmologySimulationNumberOfInitialGrids > MAX_INITIAL_GRIDS) {
    ENZO_FAIL("Too many InitialGrids! increase MAX_INITIAL_GRIDS\n");
  }
 
  if (CosmologySimulationDensityName == NULL && MultiSpecies+RadiativeCooling > 0) {
    fprintf(stderr, "warning: no density field; setting MultiSpecies/RadiativeCooling = 0\n");
    MultiSpecies = RadiativeCooling = 0;
  }

  if (CosmologySimulationParticleVelocityNames[0] != NULL &&
      !CosmologySimulationCalculatePositions) {
    ENZO_FAIL("CosmologySimulation: 1-component files only valid for use with "
	    "CosmologySimulationCalculatePositions.\n");
  }
  // If temperature is left unset, set it assuming that T=550 K at z=200
 
  if (CosmologySimulationInitialTemperature == FLOAT_UNDEFINED)
    CosmologySimulationInitialTemperature = 550.0 *
      POW((1.0 + InitialRedshift)/(1.0 + 200), 2);
 

  /* Convert from Gauss */
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1, PressureUnits=1.,MagneticUnits=1., a=1,dadt=0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  PressureUnits = DensityUnits * VelocityUnits*VelocityUnits;
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    CosmologySimulationInitialUniformBField[dim] /= MagneticUnits;
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("magnetic field: dim %"ISYM", %"FSYM" %"ESYM" \n", dim, MagneticUnits, 
	     CosmologySimulationInitialUniformBField[dim]);
  }
  // Generate the grids and set-up the hierarchy
 
  HierarchyEntry *GridsList[MAX_INITIAL_GRIDS];
  GridsList[0] = &TopGrid;
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Start loop creating initial grid hierarchy entries, from one\n");
 
  for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids; gridnum++) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("Create hierarchy entry for initial grid %"ISYM"\n", gridnum);
 
    // Create a spot in the hierarchy
 
    Subgrid    = new HierarchyEntry;
 
    // Find where to put this new grid
 
    int ParentGrid = INT_UNDEFINED;
 
    for (i = 0; i < gridnum; i++)
      if (CosmologySimulationGridLevel[i] ==
	  CosmologySimulationGridLevel[gridnum]-1)
	for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	  if (CosmologySimulationGridLeftEdge[gridnum][dim] <
	      CosmologySimulationGridLeftEdge[i][dim]       ||
	      CosmologySimulationGridRightEdge[gridnum][dim] >
	      CosmologySimulationGridRightEdge[i][dim]       )
	    break;
	  ParentGrid = i;
	}
 
    if (ParentGrid == INT_UNDEFINED) {
      ENZO_VFAIL("Grid %"ISYM" has no valid parent.\n", gridnum)
    }
 
    // Insert this grid at the appropriate position in the subgrid chain
 
    GridsList[gridnum] = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = GridsList[ParentGrid]->NextGridNextLevel;
    Subgrid->ParentGrid        = GridsList[ParentGrid];
    GridsList[ParentGrid]->NextGridNextLevel = Subgrid;
 
    // Error check for consistency and add ghost zones to dimension
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      FLOAT SubgridCellSize = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
	FLOAT(MetaData.TopGridDims[dim]*
	      POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));

      if (debug1) {
	printf("  %"GSYM"\n", SubgridCellSize);
	printf("  %"GSYM"\n", POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));
	printf("  %"ISYM" %"ISYM"\n", MetaData.TopGridDims[dim], dim);
	printf("  %"GSYM" %"GSYM"\n", CosmologySimulationGridRightEdge[gridnum][dim],
	       CosmologySimulationGridLeftEdge[gridnum][dim]);
	printf("  %"ISYM"\n", nint((CosmologySimulationGridRightEdge[gridnum][dim] -
				    CosmologySimulationGridLeftEdge[gridnum][dim]   )
				   /SubgridCellSize));
      } // ENDIF debug1
 
      // Check if declared size matches left/right edges
 
      if (nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[gridnum][dim]   )
	       /SubgridCellSize) != CosmologySimulationGridDimension[gridnum][dim]) {
	fprintf(stderr, "Subgrid inconsistency: grid %"ISYM", dim %"ISYM"\n",
		gridnum, dim);
	fprintf(stderr, " subgrid: %"GOUTSYM" -> %"GOUTSYM", CellSize = %"GOUTSYM"\n",
	      CosmologySimulationGridLeftEdge[gridnum][dim],
	      CosmologySimulationGridRightEdge[gridnum][dim], SubgridCellSize);
	ENZO_FAIL("Subgrid Inconsistency!\n");
      }
 
      // Check if left/right edge fall on Parent cell boundary
 
      if (nint((CosmologySimulationGridLeftEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ||
	  nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ) {
	fprintf(stderr, "Subgrid inconsistency: grid %"ISYM", dim %"ISYM"\n",
		gridnum, dim);
	fprintf(stderr, "left or right edges are not on parent cell edge.\n");
	ENZO_FAIL("Subgrid Inconsistency!\n");
      }
 
      // Add ghost zones
 
      CosmologySimulationGridDimension[gridnum][dim] += 2*DEFAULT_GHOST_ZONES;
 
    } // end of loop over dimensions
 
 
    // Create a new subgrid and initialize it
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(Subgrid->ParentGrid->GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank,
				   CosmologySimulationGridDimension[gridnum],
				   CosmologySimulationGridLeftEdge[gridnum],
				   CosmologySimulationGridRightEdge[gridnum],
				   0);
 
    // If subgrids are static, convert to static regions
 
    if (CosmologySimulationSubgridsAreStatic == TRUE) {
      for (region = 0; region < MAX_STATIC_REGIONS; region++)
	if (StaticRefineRegionLevel[region] == INT_UNDEFINED) {
	  StaticRefineRegionLevel[region] =
	    CosmologySimulationGridLevel[gridnum] - 1;
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] =
	      CosmologySimulationGridLeftEdge[gridnum][dim];
	    StaticRefineRegionRightEdge[region][dim] =
	      CosmologySimulationGridRightEdge[gridnum][dim];
	  }
	  for (dim = MetaData.TopGridRank; dim < MAX_DIMENSION; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] = DomainLeftEdge[dim];
	    StaticRefineRegionRightEdge[region][dim] = DomainRightEdge[dim];
	  }
	  break;
	}
      if (region == MAX_STATIC_REGIONS) {
	ENZO_FAIL("Increase number of static refine regions\n");
      }
    }
 
    // Remove ghost zones from dim
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      CosmologySimulationGridDimension[gridnum][dim] -= 2*DEFAULT_GHOST_ZONES;
 
  } // end: loop over gridnums
 
 
  //---------------------------------------------------------------------------
 
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Start loop initializing generated grids, from zero\n");
 
  // Initialize the previously-generated grids
 
  for (gridnum = 0; gridnum < CosmologySimulationNumberOfInitialGrids; gridnum++) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("RH: CosmologySimulation: Initializing grid %"ISYM"\n", gridnum);
 
    // If there is more than one grid, add the grid number to the name
 
    if (CosmologySimulationNumberOfInitialGrids > 1) {
 
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("CosmologySimulation: Initializing grid %"ISYM"\n", gridnum);
 
      if (CosmologySimulationDensityName)
	sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationDensityName, gridnum);
      if (CosmologySimulationTotalEnergyName)
	sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationTotalEnergyName, gridnum);
      if (CosmologySimulationGasEnergyName)
	sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationGasEnergyName, gridnum);
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (CosmologySimulationVelocityNames[dim])
	  sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationVelocityNames[dim], gridnum);
      if (CosmologySimulationParticlePositionName)
	sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticlePositionName, gridnum);
      if (CosmologySimulationParticleVelocityName)
	sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleVelocityName, gridnum);
      if (CosmologySimulationParticleDisplacementName)
	sprintf(ParticleDisplacementName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleDisplacementName, gridnum);
      if (CosmologySimulationParticleMassName)
	sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleMassName, gridnum);
      if (CosmologySimulationParticleTypeName)
        sprintf(ParticleTypeName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
                CosmologySimulationParticleTypeName, gridnum);
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	if (CosmologySimulationParticlePositionNames[dim])
	  sprintf(ParticlePositionNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationParticlePositionNames[dim], gridnum);
	if (CosmologySimulationParticleVelocityNames[dim])
	  sprintf(ParticleVelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationParticleVelocityNames[dim], gridnum);
	if (CosmologySimulationParticleDisplacementNames[dim])
	  sprintf(ParticleDisplacementNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationParticleDisplacementNames[dim], gridnum);
      }
	  
 
    } else {
 
      DensityName            = CosmologySimulationDensityName;
      TotalEnergyName        = CosmologySimulationTotalEnergyName;
      GasEnergyName          = CosmologySimulationGasEnergyName;
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
      ParticlePositionName   = CosmologySimulationParticlePositionName;
      ParticleVelocityName   = CosmologySimulationParticleVelocityName;
      ParticleDisplacementName = CosmologySimulationParticleDisplacementName;
      ParticleMassName       = CosmologySimulationParticleMassName;
      ParticleTypeName       = CosmologySimulationParticleTypeName;
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	ParticlePositionNames[dim] = CosmologySimulationParticlePositionNames[dim];
	ParticleVelocityNames[dim] = CosmologySimulationParticleVelocityNames[dim];
	ParticleDisplacementNames[dim] = CosmologySimulationParticleDisplacementNames[dim];
      }

    }
 
    // If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
    // otherwise just set to false
 
    SubgridsAreStatic = (GridsList[gridnum]->NextGridNextLevel == NULL) ?
      FALSE : CosmologySimulationSubgridsAreStatic;
 
    // Initialize the grid by reading in (no) data
 
    int TotalRefinement = nint(POW(FLOAT(RefineBy),
				   CosmologySimulationGridLevel[gridnum]));
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("Call CSIG for gridnum %"ISYM" with TR %"ISYM" and Dname %s\n", gridnum, TotalRefinement, DensityName);
 
    if (GridsList[gridnum]->GridData->NestedCosmologySimulationInitializeGrid(
			     gridnum,
			     CosmologySimulationOmegaBaryonNow,
			       CosmologySimulationOmegaCDMNow,
			       CosmologySimulationInitialTemperature,
			     DensityName, TotalEnergyName,
			       GasEnergyName, VelocityNames,
			       ParticlePositionName, ParticleVelocityName, 
			       ParticleDisplacementName,
			       ParticleMassName, ParticleTypeName, 
			       ParticleVelocityNames, ParticleDisplacementNames,
			     SubgridsAreStatic, TotalRefinement,
			     CosmologySimulationInitialFractionHII,
			     CosmologySimulationInitialFractionHeII,
			     CosmologySimulationInitialFractionHeIII,
			     CosmologySimulationInitialFractionHM,
			     CosmologySimulationInitialFractionH2I,
			     CosmologySimulationInitialFractionH2II,
			     CosmologySimulationInitialFractionMetal,
			     CosmologySimulationUseMetallicityField,
			     MetaData.NumberOfParticles,
			     CosmologySimulationManuallySetParticleMassRatio,
			     CosmologySimulationManualParticleMassRatio,
			     CosmologySimulationCalculatePositions,
			     CosmologySimulationGridLeftEdge[gridnum],
			     CosmologySimulationGridRightEdge[gridnum],
			     CosmologySimulationInitialUniformBField
						       ) == FAIL) {
      ENZO_FAIL("Error in grid->NestedCosmologySimulationInitializeGrid.\n");
    }
 
    // Set boundary conditions if necessary
 
  } // end loop over initial grids
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("End of loop over initial grids\n");
 
  //---------------------------------------------------------------------------
 
 
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set).
     Note: multiply MinimumMassForRefinement by the OmegaBaryonNow since the
     routine that uses this parameter only counts baryonic mass. */
 
  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
 
      MinimumMassForRefinement[i] = CosmologySimulationOmegaBaryonNow/
	                            OmegaMatterNow;
      if (CellFlaggingMethod[i] == 4)
	MinimumMassForRefinement[i] = CosmologySimulationOmegaCDMNow/
	                              OmegaMatterNow;
 
      MinimumMassForRefinement[i] *= MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *=
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }
  }
 
  // set up field names and units
 
  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = Vel1Name;
  if (MetaData.TopGridRank > 1 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    DataLabel[i++] = Vel2Name;
  if (MetaData.TopGridRank > 2 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    DataLabel[i++] = Vel3Name;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
    if(UseDivergenceCleaning){
      DataLabel[i++] = Phi_pName;
      DataLabel[i++] = DebugName;
    }
  }
   if (MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }
 
  if (CosmologySimulationUseMetallicityField) {
    DataLabel[i++] = MetalName;
    if(MultiMetals){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
  }
 
  if(STARMAKE_METHOD(COLORED_POP3_STAR)){
    DataLabel[i++] = ForbidName;
  }

  if (WritePotential)
    DataLabel[i++] = GPotName;
 
  if (CRModel) {
    DataLabel[i++] = MachName;
    if(StorePreShockFields){
      DataLabel[i++] = PSTempName;
      DataLabel[i++] = PSDenName;
    }
    DataLabel[i++] = CRName;
  } 
 

  for (j = 0; j < i; j++)
    DataUnits[j] = NULL;
 
  // Write parameters to parameter output file
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CosmologySimulationOmegaBaryonNow       = %"FSYM"\n",
	    CosmologySimulationOmegaBaryonNow);
    fprintf(Outfptr, "CosmologySimulationOmegaCDMNow          = %"FSYM"\n",
	    CosmologySimulationOmegaCDMNow);
    fprintf(Outfptr, "CosmologySimulationInitialTemperature   = %"FSYM"\n\n",
	    CosmologySimulationInitialTemperature);
 
    fprintf(Outfptr, "CosmologySimulationDensityName          = %s\n",
	    CosmologySimulationDensityName);
    if (CosmologySimulationTotalEnergyName)
    fprintf(Outfptr, "CosmologySimulationTotalEnergyName      = %s\n",
	    CosmologySimulationTotalEnergyName);
    if (CosmologySimulationGasEnergyName)
    fprintf(Outfptr, "CosmologySimulationGasEnergyName        = %s\n",
	    CosmologySimulationGasEnergyName);
    fprintf(Outfptr, "CosmologySimulationVelocity1Name        = %s\n",
	    CosmologySimulationVelocityNames[0]);
    fprintf(Outfptr, "CosmologySimulationVelocity2Name        = %s\n",
	    CosmologySimulationVelocityNames[1]);
    fprintf(Outfptr, "CosmologySimulationVelocity3Name        = %s\n",
	    CosmologySimulationVelocityNames[2]);
    if (CosmologySimulationParticlePositionName)
    fprintf(Outfptr, "CosmologySimulationParticlePositionName = %s\n",
	    CosmologySimulationParticlePositionName);
    if (CosmologySimulationParticleVelocityName)
      fprintf(Outfptr, "CosmologySimulationParticleVelocityName = %s\n",
	      CosmologySimulationParticleVelocityName);
    if (CosmologySimulationParticleDisplacementName)
	fprintf(Outfptr, "CosmologySimulationParticleDisplacementName = %s\n",
	    CosmologySimulationParticleDisplacementName);
    if (CosmologySimulationParticleMassName)
      fprintf(Outfptr, "CosmologySimulationParticleMassName     = %s\n",
	      CosmologySimulationParticleMassName);
    if (CosmologySimulationParticleTypeName)
      fprintf(Outfptr, "CosmologySimulationParticleTypeName     = %s\n\n",
	      CosmologySimulationParticleTypeName);
    fprintf(Outfptr, "CosmologySimulationParticleVelocity1Name = %s\n",
	    CosmologySimulationVelocityNames[0]);
    fprintf(Outfptr, "CosmologySimulationParticleVelocity2Name = %s\n",
	    CosmologySimulationVelocityNames[1]);
    fprintf(Outfptr, "CosmologySimulationParticleVelocity3Name = %s\n",
	    CosmologySimulationVelocityNames[2]);
    if (CosmologySimulationParticleDisplacementNames) {	
      fprintf(Outfptr, "CosmologySimulationParticleDisplacement1Name = %s\n",
	      CosmologySimulationParticleDisplacementNames[0]);
      fprintf(Outfptr, "CosmologySimulationParticleDisplacement2Name = %s\n",
	      CosmologySimulationParticleDisplacementNames[1]);
      fprintf(Outfptr, "CosmologySimulationParticleDisplacement3Name = %s\n",
	      CosmologySimulationParticleDisplacementNames[2]);
    }
    fprintf(Outfptr, "CosmologySimulationNumberOfInitialGrids = %"ISYM"\n",
	    CosmologySimulationNumberOfInitialGrids);
    fprintf(Outfptr, "CosmologySimulationSubgridsAreStatic    = %"ISYM"\n",
	    CosmologySimulationSubgridsAreStatic);
    fprintf(Outfptr, "CosmologySimulationCalculatePositions   = %"ISYM"\n",
	    CosmologySimulationCalculatePositions);
    
 
    for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids;
	 gridnum++) {
      fprintf(Outfptr, "CosmologySimulationGridLeftEdge[%"ISYM"]     = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridLeftEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridRightEdge[%"ISYM"]    = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridRightEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridDimension[%"ISYM"]    = ", gridnum);
      WriteListOfInts(Outfptr, MetaData.TopGridRank,
		      CosmologySimulationGridDimension[gridnum]);
    }
 
    fprintf(Outfptr, "\n");
 
    fprintf(Outfptr, "CosmologySimulationInitialFractionHII   = %"GSYM"\n",
	    CosmologySimulationInitialFractionHII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeII  = %"GSYM"\n",
	    CosmologySimulationInitialFractionHeII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeIII = %"GSYM"\n",
	    CosmologySimulationInitialFractionHeIII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHM    = %"GSYM"\n",
	    CosmologySimulationInitialFractionHM);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2I   = %"GSYM"\n",
	    CosmologySimulationInitialFractionH2I);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2II  = %"GSYM"\n",
	    CosmologySimulationInitialFractionH2II);
    fprintf(Outfptr, "CosmologySimulationInitialFractionMetal = %"GSYM"\n",
	    CosmologySimulationInitialFractionMetal);
    fprintf(Outfptr, "CosmologySimulationUseMetallicityField  = %"ISYM"\n\n",
	    CosmologySimulationUseMetallicityField);

    float CSBField[MAX_DIMENSION];  // in proper Gauss
    for (int dim = 0; dim < MAX_DIMENSION; dim++) 
      CSBField[dim] = CosmologySimulationInitialUniformBField[dim] * MagneticUnits;
    fprintf(Outfptr, "CosmologySimulationInitialUniformBField = ");
    WriteListOfFloats(Outfptr, 3, CSBField);

  }
 
  // Clean up
 
  delete [] dummy;
 
  return SUCCESS;
}
 
 
 
 
void NestedRecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count);
 
 
// Re-call the initializer on level zero grids.
// Used in case of ParallelRootGridIO.
 
 
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
				          TopGridData &MetaData)
{
 
 
  int dim, gridnum = 0;
 
  HierarchyEntry *CurrentGrid;
  HierarchyEntry *Temp;
 
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
    *ParticlePositionName = NULL, *ParticleVelocityName = NULL,
    *ParticleDisplacementName = NULL,
    *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION],
    *ParticleTypeName = NULL, *ParticleVelocityNames[MAX_DIMENSION],
    *ParticleDisplacementNames[MAX_DIMENSION], *ParticlePositionNames[MAX_DIMENSION];
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ParticleVelocityNames[dim] = NULL;
    ParticlePositionNames[dim] = NULL;
    ParticleDisplacementNames[dim] = NULL;
    VelocityNames[dim] = NULL;
  }
 
  CurrentGrid = TopGrid;

  /* Loop over initial grids and reinitialize each one. */

  PINT ParticleCount = 0, ParticleTempCount;
  for (gridnum = 0; gridnum < CosmologySimulationNumberOfInitialGrids; gridnum++) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("NestedCosmologySimulation: ReInitializing grid %"ISYM"\n", gridnum);
 
    // If there is more than one grid, add the grid number to the name
 
    if (CosmologySimulationNumberOfInitialGrids > 1) {
 
      if (CosmologySimulationDensityName)
	sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationDensityName, gridnum);
      if (CosmologySimulationTotalEnergyName)
	sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationTotalEnergyName, gridnum);
      if (CosmologySimulationGasEnergyName)
	sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationGasEnergyName, gridnum);
 
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	if (CosmologySimulationVelocityNames[dim])
	  sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationVelocityNames[dim], gridnum);
 
      if (CosmologySimulationParticlePositionName)
	sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticlePositionName, gridnum);
      if (CosmologySimulationParticleDisplacementName)
	sprintf(ParticleDisplacementName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticleDisplacementName, gridnum);
      if (CosmologySimulationParticleVelocityName)
	sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticleVelocityName, gridnum);
      if (CosmologySimulationParticleMassName)
	sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticleMassName, gridnum);
      if (CosmologySimulationParticleTypeName)
	sprintf(ParticleTypeName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
              CosmologySimulationParticleTypeName, gridnum);
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	if (CosmologySimulationParticleVelocityNames[dim])
	  sprintf(ParticleVelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleVelocityNames[dim], gridnum);
	if (CosmologySimulationParticlePositionNames[dim])
	  sprintf(ParticlePositionNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticlePositionNames[dim], gridnum);
	if (CosmologySimulationParticleDisplacementNames[dim])
	  sprintf(ParticleDisplacementNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleDisplacementNames[dim], gridnum);
      }
 
    } else {
 
      DensityName            = CosmologySimulationDensityName;
      TotalEnergyName        = CosmologySimulationTotalEnergyName;
      GasEnergyName          = CosmologySimulationGasEnergyName;
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
      ParticlePositionName   = CosmologySimulationParticlePositionName;
      ParticleDisplacementName = CosmologySimulationParticleDisplacementName;
      ParticleVelocityName   = CosmologySimulationParticleVelocityName;
      ParticleMassName       = CosmologySimulationParticleMassName;
      ParticleTypeName       = CosmologySimulationParticleTypeName;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	ParticleVelocityNames[dim]   = CosmologySimulationParticleVelocityNames[dim];
	ParticlePositionNames[dim]   = CosmologySimulationParticlePositionNames[dim];
	ParticleDisplacementNames[dim] = CosmologySimulationParticleDisplacementNames[dim];
      }
 
    }
 
    // If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
    // otherwise just set to false
 
    int SubgridsAreStatic = (CurrentGrid->NextGridNextLevel == NULL) ?
      FALSE : CosmologySimulationSubgridsAreStatic;
 
    // Call grid initializer.  Use TotalRefinement = -1 to flag real read
 
    int TotalRefinement = -1;
 
    // Loop over all grids on this level
 
    Temp = CurrentGrid;
    while (Temp != NULL) {
      ParticleTempCount = ParticleCount; // set particle count to beginning of this level (not used for ring IO)
      if (Temp->GridData->NestedCosmologySimulationInitializeGrid
	  (gridnum, CosmologySimulationOmegaBaryonNow,
	   CosmologySimulationOmegaCDMNow,
	   CosmologySimulationInitialTemperature,
	   DensityName, TotalEnergyName,
	   GasEnergyName, VelocityNames,
	   ParticlePositionName, ParticleVelocityName, 
	   ParticleDisplacementName,
	   ParticleMassName, ParticleTypeName, 
	   ParticleVelocityNames, ParticleDisplacementNames,
	   SubgridsAreStatic, TotalRefinement,
	   CosmologySimulationInitialFractionHII,
	   CosmologySimulationInitialFractionHeII,
	   CosmologySimulationInitialFractionHeIII,
	   CosmologySimulationInitialFractionHM,
	   CosmologySimulationInitialFractionH2I,
	   CosmologySimulationInitialFractionH2II,
	   CosmologySimulationInitialFractionMetal,
	   CosmologySimulationUseMetallicityField,
	   ParticleTempCount,
	   CosmologySimulationManuallySetParticleMassRatio,
	   CosmologySimulationManualParticleMassRatio,
	   CosmologySimulationCalculatePositions,
	   CosmologySimulationGridLeftEdge[gridnum],
	   CosmologySimulationGridRightEdge[gridnum],
	   CosmologySimulationInitialUniformBField
	   ) == FAIL) {
	ENZO_FAIL("Error in grid->NestedCosmologySimulationInitializeGrid.\n");
      }
 
      Temp = Temp->NextGridThisLevel;
    } // end: loop over grids on this level

    /* Once we have read in all the grids, update the current particle count
       by the number returned, which is the total number of particles on that level
       (this only works for one initial sub-region per level). */

    ParticleCount = ParticleTempCount; // set particle count (not used for ring IO)

    // Go down to the grid(s) on the next level
 
    CurrentGrid = CurrentGrid->NextGridNextLevel;
 
  } // end loop over initial grid levels

  // Create tracer particles (on top grid)
  
  if (TracerParticleCreationSpacing > 0) {

    // Not sure if this works.
    if (ParallelRootGridIO) 
      printf("Warning: Tracer particles with parallel IO not tested.\n");
 
    PINT DummyNumberOfParticles = 0;
 
    Temp = TopGrid;
    while (Temp != NULL) {
      if (Temp->GridData->TracerParticleCreateParticles(
                                     TracerParticleCreationLeftEdge,
                                     TracerParticleCreationRightEdge,
                                     TracerParticleCreationSpacing,
                                     DummyNumberOfParticles) == FAIL) {
	ENZO_FAIL("Error in grid->TracerParticleCreateParticles\n");
      }
 
      Temp = Temp->NextGridThisLevel;
    }

  } // end: if (TracerParticleCreationSpacing > 0)
 
  // Get the global particle count
 
  int LocalNumberOfParticles;;
  ParticleCount = 0;
 
  CurrentGrid = TopGrid;
  while (CurrentGrid != NULL) {
 
    Temp = CurrentGrid;
    while (Temp != NULL) {
 
      LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
      // printf("OldLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles );
 
#ifdef USE_MPI
      CommunicationAllReduceValues(&LocalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */
      Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);
 
      //LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
      // printf("NewLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles );
      ParticleCount += LocalNumberOfParticles;
 
      Temp = Temp->NextGridThisLevel;
    }
 
    CurrentGrid = CurrentGrid->NextGridNextLevel;
 
  }
 
  // Loop over grids and set particle ID number.
  // This is done for ring IO but is not necessary for regular IO
  //  (particle IDs are set during read and do not depend on grid distribution)

  if (ParallelParticleIO || CosmologySimulationCalculatePositions) {
    Temp = TopGrid;
    ParticleCount = 0;
 
    NestedRecursivelySetParticleCount(Temp, &ParticleCount);

  }
 
  if (debug)
    printf("FinalParticleCount = %"PISYM"\n", ParticleCount);
 
  MetaData.NumberOfParticles = ParticleCount;

  /* Now set the grids' parents correctly instead of having all of the
     siblings have the same parent.  Before this left the work for
     RebuildHierarchy, which will send all of the regions and
     particles back to one processor! */

  bool match;
  int level, Rank, TempDims[MAX_DIMENSION];
  FLOAT GridCenter[MAX_DIMENSION];
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  LevelHierarchyEntry *CArray, *PArray;
  HierarchyEntry *Current, *Parent;  // For convenience
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT LeftParent[MAX_DIMENSION], RightParent[MAX_DIMENSION];

  // Level arrays are convenient
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  AddLevel(LevelArray, TopGrid, 0);

  /* Now that we have the pointers stored in LevelArray, let's reset
     all of the hierarchy pointers (except for the level-0
     NextGridThisLevel) */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (CArray = LevelArray[level]; CArray; CArray = CArray->NextGridThisLevel) {
      if (level > 0)
	CArray->GridHierarchyEntry->NextGridThisLevel = NULL;
      CArray->GridHierarchyEntry->NextGridNextLevel = NULL;
      CArray->GridHierarchyEntry->ParentGrid = NULL;
    } // ENDFOR grids
  } // ENDFOR level

  // Loop over all child grids, looking for parents
  for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (CArray = LevelArray[level]; CArray; CArray = CArray->NextGridThisLevel) {

      Current = CArray->GridHierarchyEntry;
      Current->GridData->ReturnGridInfo(&Rank, TempDims, LeftEdge, RightEdge);
      for (dim = 0; dim < Rank; dim++)
	GridCenter[dim] = 0.5*(LeftEdge[dim] + RightEdge[dim]);

      // Look for parents in level-1
      for (PArray = LevelArray[level-1]; PArray; 
	   PArray = PArray->NextGridThisLevel) {

	Parent = PArray->GridHierarchyEntry;
	Parent->GridData->ReturnGridInfo(&Rank, TempDims, LeftParent, RightParent);

	// Grid won't necessarily be fully contained in initial
	// "parent".  Here we just look for overlap if the partitioned
	// grids weren't created to be fully contained within a
	// partition of the level-1 grid.

	match = true;
	for (dim = 0; dim < Rank; dim++)
	  match &= (GridCenter[dim] > LeftParent[dim]) &&
	    (GridCenter[dim] < RightParent[dim]);

	if (match) {
	  Current->ParentGrid = Parent;
	  Current->NextGridThisLevel = Parent->NextGridNextLevel;
	  Parent->NextGridNextLevel = Current;
	  break;
	} // ENDIF match

      } // ENDFOR parent level grids
    } // ENDFOR current level grids
  } // ENDFOR level

  return SUCCESS;
}
 
 
 
 
void NestedRecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count)
{
  // Add Count to the particle id's on this grid (which start from zero
  // since we are doing a parallel root grid i/o)
 
  GridPoint->GridData->AddToParticleNumber(Count);
 
  // Recursively apply this to siblings and children
 
  if (GridPoint->NextGridThisLevel != NULL)
    NestedRecursivelySetParticleCount(GridPoint->NextGridThisLevel, Count);
 
  if (GridPoint->NextGridNextLevel != NULL)

    NestedRecursivelySetParticleCount(GridPoint->NextGridNextLevel, Count);
 
  CommunicationBroadcastValue(Count, ROOT_PROCESSOR);
  return;
}
