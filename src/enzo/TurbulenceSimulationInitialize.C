/***********************************************************************
/
/  INITIALIZE A TURBULENCE SIMULATION
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:  Initialize a turbulence simulation.  Reads in initial data
/            for the root grid.
/
/  DEFAULT:  Quasi-isothermal forced turbulence.
/            Requires RAREFACTION1 (or 0), defined in euler.src.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

 
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
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
#include "fortran.def"
#include "error.h"
#include "message.h"

const char config_turbulence_simulation_defaults[] =
"### TURBULENCE SIMULATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    TurbulenceSimulation: {\n"
"        InitialDensity				= -99999.0;\n"
"        InitialDensityPerturbationAmplitude 	= -99999.0;\n"
"        InitialTemperature			= -99999.0;\n"
"        InitialPressure			= -99999.0;\n"
"        InitialMagneticField[3] 		= [5.0,5.0,5.0];\n"
"\n"
"        DensityName		= \"\";\n"
"        TotalEnergyName	= \"\";\n"
"        GasPressureName	= \"\";\n"
"        GasEnergyName		= \"\";\n"
"        VelocityNames		= \"\";\n"
"        RandomForcingNames	= \"\";\n"
"        MagneticNames		= \"\";\n"
"        SubgridsAreStatic	= TRUE;
"        NumberofInitialSubgrids= 1;
"    };\n"
"};\n";

 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationAllSumIntegerValues(int *Values, int Number);
 
/* Turbulence Parameters (that need to be shared) */
 
static float TurbulenceSimulationInitialDensityPerturbationAmplitude;
static float TurbulenceSimulationInitialTemperature;
static float TurbulenceSimulationInitialPressure;
static float TurbulenceSimulationInitialMagneticField[3];
 
static char *TurbulenceSimulationDensityName;
static char *TurbulenceSimulationTotalEnergyName;
static char *TurbulenceSimulationGasPressureName;
static char *TurbulenceSimulationGasEnergyName;
static char *TurbulenceSimulationVelocityNames[MAX_DIMENSION];
static char *TurbulenceSimulationRandomForcingNames[MAX_DIMENSION];
static char *TurbulenceSimulationMagneticNames[MAX_DIMENSION];
static int   TurbulenceSimulationSubgridsAreStatic;
static int   TurbulenceSimulationNumberOfInitialGrids;
 
 
#define MAX_INITIAL_GRIDS 10
 
int TurbulenceSimulationInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *Drive1Name = "DrivingField1";
  char *Drive2Name = "DrivingField2";
  char *Drive3Name = "DrivingField3";
  char *GravPotName = "GravPotential";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
 
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasPressureName = NULL, *GasEnergyName = NULL,
       *VelocityNames[MAX_DIMENSION],
       *RandomForcingNames[MAX_DIMENSION],
       *MagneticNames[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    RandomForcingNames[dim] = NULL;
    MagneticNames[dim] = NULL;
  }
 
  /* Set default parameters and names */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    TurbulenceSimulationVelocityNames[dim]       = NULL;
 
  int   TurbulenceSimulationGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   TurbulenceSimulationGridLevel[MAX_INITIAL_GRIDS];
  FLOAT TurbulenceSimulationGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT TurbulenceSimulationGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    TurbulenceSimulationGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    TurbulenceSimulationGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    TurbulenceSimulationGridRightEdge[0][dim] = DomainRightEdge[dim];
    TurbulenceSimulationGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }

  int InitialMagneticFieldDefined = FALSE;
  TurbulenceSimulationGridLevel[0] = 0;
 
  float TurbulenceSimulationSoundSpeed = 1.0;

  /* Error check. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (DualEnergyFormalism == TRUE && HydroMethod != Zeus_Hydro)
      fprintf(stderr, "TurbulenceSimulation: DualEnergyFormalism is ON.\n");
    if (!SelfGravity)
      fprintf(stderr, "TurbulenceSimulation: SelfGravity is OFF\n");
    if (RandomForcing && HydroMethod == Zeus_Hydro) {
      fprintf(stderr, "RandomForcing for Zeus hydro is not tested!\n");
      // ERROR_MESSAGE;
    }
  }
 
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  char dummy_arr[MAX_DIMENSION][MAX_LINE_LENGTH];


  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_turbulence_simulation_defaults);

 
    /* Read parameters */

     
    Param.GetScalar(dummy, "Problem.TurbulenceSimulation.DensityName");
    if( strlen(dummy) > 0 ) {
      TurbulenceSimulationDensityName = new char[MAX_LINE_LENGTH];
      strcpy(TurbulenceSimulationDensityName,dummy);
    }

    Param.GetScalar(dummy, "Problem.TurbulenceSimulation.TotalEnergyName");
    if( strlen(dummy) > 0 ) {
      TurbulenceSimulationTotalEnergyName = new char[MAX_LINE_LENGTH];
      strcpy(TurbulenceSimulationTotalEnergyName,dummy);
    }

    Param.GetScalar(dummy, "Problem.TurbulenceSimulation.GasPressureName");
    if( strlen(dummy) > 0 ) {
      TurbulenceSimulationGasPressureName = new char[MAX_LINE_LENGTH];
      strcpy(TurbulenceSimulationGasPressureName,dummy);
    }

    Param.GetScalar(dummy, "Problem.TurbulenceSimulation.GasEnergyName");
    if( strlen(dummy) > 0 ) {
      TurbulenceSimulationGasEnergyName = new char[MAX_LINE_LENGTH];
      strcpy(TurbulenceSimulationGasEnergyName,dummy);
    }

    Param.GetArray(dummy_arr, "Problem.TurbulenceSimulation.VelocityNames");
    for( i=0; i<MAX_DIMENSION; i++) {
      if( strlen(dummy_arr[i]) > 0 ) {
        TurbulenceSimulationVelocityNames[i] = new char[MAX_LINE_LENGTH];
        strcpy(TurbulenceSimulationVelocityNames[i], dummy_arr[i]);
      }
    }

    Param.GetArray(dummy_arr, "Problem.TurbulenceSimulation.RandomForcingNames");
    for( i=0; i<MAX_DIMENSION; i++) {
      if( strlen(dummy_arr[i]) > 0 ) {
        TurbulenceSimulationRandomForcingNames[i] = new char[MAX_LINE_LENGTH];
        strcpy(TurbulenceSimulationRandomForcingNames[i], dummy_arr[i]);
      }
    }

    Param.GetArray(dummy_arr, "Problem.TurbulenceSimulation.MagneticNames");
    for( i=0; i<MAX_DIMENSION; i++) {
      if( strlen(dummy_arr[i]) > 0 ) {
        TurbulenceSimulationMagneticNames[i] = new char[MAX_LINE_LENGTH];
        strcpy(TurbulenceSimulationMagneticNames[i], dummy_arr[i]);
      }
    }

  Param.GetScalar(TurbulenceSimulationInitialTemperature, "Problem.TurbulenceSimulationInitialTemperature");
  Param.GetScalar(TurbulenceSimulationInitialDensity, "Problem.TurbulenceSimulationInitialDensity");
  Param.GetScalar(TurbulenceSimulationSoundSpeed, "Problem.TurbulenceSimulationSoundSpeed");
  Param.GetScalar(TurbulenceSimulationInitialPressure, "Problem.TurbulenceSimulationInitialPressure");
  Param.GetScalar(TurbulenceSimulationInitialDensityPerturbationAmplitude, 
		"Problem.TurbulenceSimulation.InitialDensityPerturbationAmplitude");
  Param.GetScalar(TurbulenceSimulationNumberOfInitialGrids, 
		"Problem.TurbulenceSimulation.NumberOfInitialGrids");
  Param.GetScalar(TurbulenceSimulationSubgridsAreStatic,
		"Problem.TurbulenceSimulation.SubgridsAreStatic");


  int NumberOfInitialSubgrids = Param.Size("Problem.TurbulenceSimulation.InitialSubgrids");
  if (NumberOfInitialSubgrids+1 > MAX_INITIAL_GRIDS) {
    ENZO_VFAIL("You've exceeded the maximum number of TurbulenceSimulation initial grids (%d)!\n",
      MAX_INITIAL_GRIDS)
  }// end if

  TurbulenceSimulationNumberOfInitialGrids = NumberOfInitialSubgrids + 1;  // add 1 for the root grid

  char InitialSubgridNames[MAX_LINE_LENGTH][MAX_INITIAL_GRIDS];
  Param.GetArray(InitialSubgridNames,"Problem.TurbulenceSimulation.InitialSubgrids");

  for (i = 0; i < NumberOfInitialSubgrids; i++) {
    Param.GetArray(TurbulenceSimulationGridLeftEdge[i+1], 
        "Problem.TurbulenceSimulation.%s.LeftEdge",InitialSubgridNames[i]);
    Param.GetArray(TurbulenceSimulationGridRightEdge[i+1], 
        "Problem.TurbulenceSimulation.%s.RightEdge",InitialSubgridNames[i]);
    Param.GetArray(TurbulenceSimulationGridDimension[i+1],
        "Problem.TurbulenceSimulation.%s.Dimension",InitialSubgridNames[i]);
    Param.GetScalar(TurbulenceSimulationGridLevel[i+1], 
        "Problem.TurbulenceSimulation.%s.Level",InitialSubgridNames[i]);
  }// end i for

  Param.GetArray(TurbulenceSimulationInitialMagneticField,"Problem.TurbulenceSimulation.InitialMagneticField");
  Param.GetScalar(InitialMagneticFieldDefined,
  		"Problem.TurbulenceSimulation.MagneticFieldDefined");

    /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
  /* More error checking. */
  /* dcc removed in favor of in house generation.
  if (TurbulenceSimulationVelocityNames[0] == NULL) {
    ENZO_FAIL("Missing initial data.\n");
  }
  */
  if (CellFlaggingMethod[0] != 3)
      fprintf(stderr, "TurbulenceSimulation: check CellFlaggingMethod.\n");
 
  /* If density(temperature, actually, c^2) is left unset, set it = 1. */
 
  if (TurbulenceSimulationInitialDensity == FLOAT_UNDEFINED)
    TurbulenceSimulationInitialDensity = 1.0;
  if (TurbulenceSimulationInitialDensityPerturbationAmplitude == FLOAT_UNDEFINED)
    TurbulenceSimulationInitialDensityPerturbationAmplitude= 0.0;
  if (TurbulenceSimulationInitialTemperature == FLOAT_UNDEFINED &&
		  TurbulenceSimulationInitialPressure == FLOAT_UNDEFINED){
    TurbulenceSimulationInitialTemperature = 1.0;
  }

  if( HydroMethod == MHD_RK) {
    if( InitialMagneticFieldDefined != TRUE) {
      TurbulenceSimulationInitialMagneticField[0] = 1e-8;
      TurbulenceSimulationInitialMagneticField[1] = 1e-8;
      TurbulenceSimulationInitialMagneticField[2] = 1e-8;
    }
  }
  
  /* Check/define RandomForcing parameters [Mac Low 1999, ApJ 524, 169]. */
 
  if (RandomForcing)
    if (RandomForcingMachNumber <= 0.0) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Warning: RandomForcing is OFF\n");
      RandomForcing = 0;
    }
 
  /* If RandomForcingEdot is not set in the parameter file, get it
     from [MacLow1999] formula. Note: his formula is calibrated for
     general random forcing fields; coefficient 0.81 can potentially
     be inappropriate for a purely solenoidal forcing; also our
     Gamma is not quite 1.0. */
 
  if (RandomForcing && RandomForcingEdot < 0.0) {
    float BoxSize = TurbulenceSimulationGridRightEdge[0][0] -
                    TurbulenceSimulationGridLeftEdge[0][0];
    float BoxMass = (BoxSize*BoxSize*BoxSize)*
                    TurbulenceSimulationInitialDensity;
    float Vrms    = RandomForcingMachNumber/
      sqrt(1);//Sound speed is one.
    RandomForcingEdot = 0.81/BoxSize*BoxMass*Vrms*Vrms*Vrms;
 
  /* Approximate correction to the MacLow's factor (see eqs (7) - (8))
     for **this PPM implementation**. Seems to be OK for 64^3, 128^3 and 256^3
     Mach=3,6,10 simulations of **solenoidally** driven turbulence. */
 
    RandomForcingEdot *= 0.8;
 
  }
 
  if (RandomForcing)
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("RandomForcingEdot: %"GSYM"\n", RandomForcingEdot);
 
  /* -------------------------------------------------------------------- */
  /* Generate the root grid and set-up the hierarchy. */
 
  HierarchyEntry *GridsList;
  GridsList = &TopGrid;
 
  /* Initialize the root grid. */
 
  DensityName            = TurbulenceSimulationDensityName;
  TotalEnergyName        = TurbulenceSimulationTotalEnergyName;
  GasPressureName        = TurbulenceSimulationGasPressureName;
  GasEnergyName          = TurbulenceSimulationGasEnergyName;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    VelocityNames[dim]   = TurbulenceSimulationVelocityNames[dim];
    RandomForcingNames[dim] = TurbulenceSimulationRandomForcingNames[dim];
    MagneticNames[dim] = TurbulenceSimulationMagneticNames[dim];

  }
 
  /* Initialize the root grid by reading in data. */
 
  int TotalRefinement = nint(POW(FLOAT(RefineBy),
                                   TurbulenceSimulationGridLevel[gridnum]));
  if (GridsList->GridData->TurbulenceSimulationInitializeGrid(
			   TurbulenceSimulationInitialDensity,
			   TurbulenceSimulationInitialDensityPerturbationAmplitude,
			   TurbulenceSimulationInitialTemperature,
			   TurbulenceSimulationInitialPressure,
			   TurbulenceSimulationInitialMagneticField,
			   MagneticNames,
			   DensityName, TotalEnergyName,GasPressureName,
			   GasEnergyName, VelocityNames, RandomForcingNames,
			   TurbulenceSimulationSubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->TurbulenceSimulationInitializeGrid.\n");
  }
 
  /* Set boundary conditions if necessary. */
 
  // this will be done in EvolveHierarchy()
 
  /* set up field names and units */
 
  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
  }
  if (UseDrivingField) {
    DataLabel[i++] = Drive1Name;
    DataLabel[i++] = Drive2Name;
    DataLabel[i++] = Drive3Name;
  }
  if(WritePotential)
      DataLabel[i++] = GravPotName;

  for (j = 0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file.
     ATTENTION: printf %s on sun fails if passed a NULL pointer. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TurbulenceSimulationInitialDensity   = %"FSYM"\n\n",
	    TurbulenceSimulationInitialTemperature);
    fprintf(Outfptr, "TurbulenceSimulationInitialDensityPerturbationAmplitude = %f\n\n", 
	    TurbulenceSimulationInitialDensityPerturbationAmplitude);
 
    fprintf(Outfptr, "TurbulenceSimulationInitialTemperature   = %"FSYM"\n\n",
	    TurbulenceSimulationInitialTemperature);
    fprintf(Outfptr, "TurbulenceSimulationInitialPressure   = %f\n\n", 
	    TurbulenceSimulationInitialPressure);
 
    if (TurbulenceSimulationDensityName)
    fprintf(Outfptr, "TurbulenceSimulationDensityName          = %s\n",
	    TurbulenceSimulationDensityName);
    if (TurbulenceSimulationTotalEnergyName)
    fprintf(Outfptr, "TurbulenceSimulationTotalEnergyName      = %s\n",
	    TurbulenceSimulationTotalEnergyName);
    if (TurbulenceSimulationGasPressureName)
    fprintf(Outfptr, "TurbulenceSimulationGasPressureName      = %s\n",
	    TurbulenceSimulationGasPressureName);
    if (TurbulenceSimulationGasEnergyName)
    fprintf(Outfptr, "TurbulenceSimulationGasEnergyName        = %s\n",
	    TurbulenceSimulationGasEnergyName);
    fprintf(Outfptr, "TurbulenceSimulationVelocity1Name        = %s\n",
	    TurbulenceSimulationVelocityNames[0]);
    fprintf(Outfptr, "TurbulenceSimulationVelocity2Name        = %s\n",
	    TurbulenceSimulationVelocityNames[1]);
    fprintf(Outfptr, "TurbulenceSimulationVelocity3Name        = %s\n",
	    TurbulenceSimulationVelocityNames[2]);
    if (TurbulenceSimulationRandomForcingNames[0])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing1Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[0]);
    if (TurbulenceSimulationRandomForcingNames[1])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing2Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[1]);
    if (TurbulenceSimulationRandomForcingNames[2])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing3Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[2]);
  }
 
  /* Clean up. */
 
  delete dummy;
 
  return SUCCESS;
}
 
 
/* --------------------------------------------------------------------
   Re-call the initializer on level zero grids.  Used in case of
   ParallelRootGridIO. */
 
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
				    TopGridData &MetaData)
{
 
  /* Declarations. */
 
  int dim, gridnum = 0;
  char *DensityName = NULL, *TotalEnergyName = NULL,*GasPressureName = NULL, *GasEnergyName = NULL,
       *ParticlePositionName = NULL, *ParticleVelocityName = NULL, 
       *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION], 
    *RandomForcingNames[MAX_DIMENSION],
    *MagneticNames[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    RandomForcingNames[dim] = NULL;
    MagneticNames[dim] = NULL;
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("TurbulenceSimulation: ReInitializing grid %"ISYM"\n", gridnum);
 
  /* If there is more than one grid, add the grid number to the name. */
 
  if (TurbulenceSimulationNumberOfInitialGrids > 1) {
 
    if (TurbulenceSimulationDensityName)
      sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationDensityName, gridnum);
    if (TurbulenceSimulationTotalEnergyName)
      sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationTotalEnergyName, gridnum);
    if (TurbulenceSimulationGasPressureName)
      sprintf(GasPressureName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationGasPressureName, gridnum);
    if (TurbulenceSimulationGasEnergyName)
      sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationGasEnergyName, gridnum);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (TurbulenceSimulationVelocityNames[dim])
	sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		TurbulenceSimulationVelocityNames[dim], gridnum);
      if (TurbulenceSimulationRandomForcingNames[dim])
	sprintf(RandomForcingNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		TurbulenceSimulationRandomForcingNames[dim], gridnum);
    }
 
  } else {
 
    DensityName            = TurbulenceSimulationDensityName;
    TotalEnergyName        = TurbulenceSimulationTotalEnergyName;
    GasPressureName        = TurbulenceSimulationGasPressureName;
    GasEnergyName          = TurbulenceSimulationGasEnergyName;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      VelocityNames[dim]   = TurbulenceSimulationVelocityNames[dim];
      RandomForcingNames[dim] = TurbulenceSimulationRandomForcingNames[dim];
      MagneticNames[dim] = TurbulenceSimulationMagneticNames[dim];

    }
 
  }
 
  /* Call grid initializer.  Use TotalRefinement = -1 to flag real read. */
 
  int TotalRefinement = -1;
 
  /* Loop over level zero grid. */
 
  HierarchyEntry *Temp = TopGrid;
  while (Temp != NULL) {
 
    if (Temp->GridData->TurbulenceSimulationInitializeGrid(
		        TurbulenceSimulationInitialDensity,
			TurbulenceSimulationInitialDensityPerturbationAmplitude,
		        TurbulenceSimulationInitialTemperature,
		        TurbulenceSimulationInitialPressure,
		        TurbulenceSimulationInitialMagneticField,
			MagneticNames,
		        DensityName, TotalEnergyName,GasPressureName,
		        GasEnergyName, VelocityNames, RandomForcingNames,
			TurbulenceSimulationSubgridsAreStatic,
			TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->TurbulenceSimulationInitializeGrid.\n");

    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}
 
