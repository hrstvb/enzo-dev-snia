/***********************************************************************
/
/  READ A PARAMETER FILE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       June 25, 2006
/  modified2:  Robert Harkness
/  date:       February 29th, 2008
/  modified3:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine reads the parameter file in the argument and sets parameters
//   based on it.

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "hydro_rk/EOS.h" 

//#include "parameters.h"


/* This variable is declared here and only used in Grid_ReadGrid. */
 


/* function prototypes */

void my_exit(int status); 
int CosmologyReadParameters(FLOAT *StopTime, FLOAT *InitTime);
int ReadUnits();
int InitializeCosmicRayData();
int InitializeRateData(FLOAT Time);
int InitializeEquilibriumCoolData(FLOAT Time);
int InitializeGadgetEquilibriumCoolData(FLOAT Time);
int InitializeRadiationFieldData(FLOAT Time);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int ReadEvolveRefineFile(void);
 
int CheckShearingBoundaryConsistency(TopGridData &MetaData); 
void get_uuid(char *buffer);

int ReadParameterFile(TopGridData &MetaData, float *Initialdt)
{
#ifndef CONFIG_USE_LIBCONFIG
	
  /* declarations */
  
  int i, dim, ret;
  float TempFloat;
  int comment_count = 0;
	
  char *dummy = new char[MAX_LINE_LENGTH];


  /* Initialize all MetaData char arrays to NULL */

  MetaData.ResubmitCommand        = NULL;
  MetaData.RestartDumpName        = NULL;
  MetaData.DataDumpName           = NULL;
  MetaData.HistoryDumpName        = NULL;
  MetaData.MovieDumpName          = NULL;
  MetaData.TracerParticleDumpName = NULL;
  MetaData.RedshiftDumpName       = NULL;
  MetaData.RestartDumpDir         = NULL;
  MetaData.DataDumpDir            = NULL;
  MetaData.HistoryDumpDir         = NULL;
  MetaData.MovieDumpDir           = NULL;
  MetaData.TracerParticleDumpDir  = NULL;
  MetaData.RedshiftDumpDir        = NULL;
  MetaData.LocalDir               = NULL;
  MetaData.GlobalDir              = NULL;
  MetaData.MetaDataIdentifier     = NULL;
  MetaData.SimulationUUID         = NULL;
  MetaData.RestartDatasetUUID     = NULL;
  MetaData.InitialConditionsUUID  = NULL;
  MetaData.BoundaryConditionName  = NULL;
#ifdef TRANSFER
  MetaData.RadHydroParameterFname = NULL;
#endif  

  /* read MetaData parameters */

  Param.GetScalar(MetaData.CycleNumber, "Internal.InitialCycleNumber");
  Param.GetScalar(MetaData.Time, "Internal.InitialTime");
  Param.GetScalar(MetaData.CPUTime, "Internal.InitialCPUTime");
  Param.GetScalar((*Initialdt), "Internal.Initialdt");
  
  Param.GetScalar(CheckpointRestart, "Internal.OutputLabeling.CheckpointRestart");
  Param.GetScalar(MetaData.StopTime, "SimulationControl.StopTime");
  Param.GetScalar(MetaData.StopCycle, "SimulationControl.StopCycle");
  Param.GetScalar(MetaData.StopSteps, "SimulationControl.StopSteps");
  Param.GetScalar(MetaData.StopCPUTime, "SimulationControl.StopCPUTime");
  Param.GetScalar(MetaData.ResubmitOn, "SimulationControl.ResubmitOn");
  
	//printf("Size of resubmitcommand = %d",Param.Size("SimulationControl.ResubmitCommand"));
	
  Param.GetScalar(dummy, "SimulationControl.ResubmitCommand");
  if( strlen(dummy) > 0 ) {
    MetaData.ResubmitCommand = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.ResubmitCommand, dummy);
  }

  Param.GetScalar(MetaData.MaximumTopGridTimeStep, "SimulationControl.MaximumTopGridTimeStep");
  
  Param.GetScalar(MetaData.TimeLastRestartDump, "Internal.OutputLabeling.TimeLastRestartDump");
  Param.GetScalar(MetaData.dtRestartDump, "OutputControl.RestartDump.dt");
  Param.GetScalar(MetaData.TimeLastDataDump, "Internal.OutputLabeling.TimeLastDataDump");
  Param.GetScalar(MetaData.dtDataDump, "OutputControl.DataDump.dt");
  Param.GetScalar(MetaData.TimeLastHistoryDump, "Internal.OutputLabeling.TimeLastHistoryDump");
  Param.GetScalar(MetaData.dtHistoryDump, "OutputControl.HistoryDump.dt");
  
  Param.GetScalar(TracerParticleOn, "SimulationControl.TracerParticleOn");
  Param.GetScalar(ParticleTypeInFile, "OutputControl.ParticleTypeInFile");
  Param.GetScalar(OutputParticleTypeGrouping, "OutputControl.OutputParticleTypeGrouping");
  Param.GetScalar(MetaData.TimeLastTracerParticleDump, "Internal.OutputLabeling.TimeLastTracerParticleDump"); //doc?
  Param.GetScalar(MetaData.dtTracerParticleDump, "OutputControl.TracerParticleDump.dt");//doc?
  Param.GetScalar(MetaData.TimeLastInterpolatedDataDump, "Internal.OutputLabeling.TimeLastInterpolatedDataDump");//doc?
  Param.GetScalar(MetaData.dtInterpolatedDataDump, "OutputControl.dtInterpolatedDataDump");//doc?
  
  
  Param.GetArray(MetaData.NewMovieLeftEdge, "OutputControl.MovieDump.LeftEdge");
  Param.GetArray(MetaData.NewMovieRightEdge, "OutputControl.MovieDump.RightEdge");
  
  Param.GetScalar(MetaData.CycleLastRestartDump, "Internal.OutputLabeling.CycleLastRestartDump"); //not used
  Param.GetScalar(MetaData.CycleSkipRestartDump, "OutputControl.CycleDump.SkipRestartDump"); //not used
  Param.GetScalar(MetaData.CycleLastDataDump, "Internal.OutputLabeling.CycleLastDataDump");
  Param.GetScalar(MetaData.CycleSkipDataDump, "OutputControl.CycleDump.SkipDataDump");
  Param.GetScalar(MetaData.CycleLastHistoryDump, "Internal.OutputLabeling.CycleLastHistoryDump");
  Param.GetScalar(MetaData.CycleSkipHistoryDump, "OutputControl.CycleDump.SkipHistoryDump");
  Param.GetScalar(MetaData.CycleSkipGlobalDataDump, "OutputControl.CycleDump.SkipGlobalDataDump");
  Param.GetScalar(MetaData.OutputFirstTimeAtLevel, "OutputControl.OutputTriggers.OutputFirstTimeAtLevel");
  Param.GetScalar(MetaData.StopFirstTimeAtLevel, "SimulationControl.StopFirstTimeAtLevel");
  
  
  /* Maximum density directed output */
  Param.GetScalar(OutputOnDensity, "OutputControl.OutputTriggers.OutputOnDensity"); 
  Param.GetScalar(StartDensityOutputs, "OutputControl.OutputTriggers.StartDensityOutputs");
  Param.GetScalar(CurrentDensityOutput, "OutputControl.OutputTriggers.CurrentDensityOutput");
  Param.GetScalar(IncrementDensityOutput, "OutputControl.OutputTriggers.IncrementDensityOutput");
  
  
  /* Subcycle directed output */
  Param.GetScalar(MetaData.SubcycleSkipDataDump, "OutputControl.CycleDump.SubcycleSkipDataDump");
  Param.GetScalar(MetaData.SubcycleLastDataDump, "Internal.OutputLabeling.SubcycleLastDataDump");
  Param.GetScalar(MetaData.SubcycleNumber, "Internal.SubcycleNumber");
  
  
  Param.GetScalar(FileDirectedOutput, "OutputControl.FileDirectedOutput");
  Param.GetScalar(WriteBinaryHierarchy, "OutputControl.WriteBinaryHierarchy");
  
  
  Param.GetScalar(MetaData.RestartDumpNumber, "Internal.OutputLabeling.RestartDumpNumber");
  Param.GetScalar(MetaData.DataDumpNumber, "Internal.OutputLabeling.DataDumpNumber");
  Param.GetScalar(MetaData.HistoryDumpNumber, "Internal.OutputLabeling.HistoryDumpNumber");
  Param.GetScalar(MetaData.TracerParticleDumpNumber, "Internal.OutputLabeling.TracerParticleDumpNumber");
  
  Param.GetScalar(dummy, "OutputControl.RestartDump.Name");
  if( strlen(dummy) > 0 ) {
    MetaData.RestartDumpName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RestartDumpName, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.RestartDump.Dir");
  if( strlen(dummy) > 0 ) {
    MetaData.RestartDumpDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RestartDumpDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.DataDump.Name");
  if( strlen(dummy) > 0 ) {
    MetaData.DataDumpName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.DataDumpName, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.DataDump.Dir");
  if( strlen(dummy) > 0 ) {
    MetaData.DataDumpDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.DataDumpDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.RedshiftDump.Name");
  if( strlen(dummy) > 0 ) {
    MetaData.RedshiftDumpName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RedshiftDumpName, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.RedshiftDump.Dir");
  if( strlen(dummy) > 0 ) {
    MetaData.RedshiftDumpDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RedshiftDumpDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.TracerParticleDump.Name");
  if( strlen(dummy) > 0 ) {
    MetaData.TracerParticleDumpName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.TracerParticleDumpName, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.TracerParticleDump.Dir");
  if( strlen(dummy) > 0 ) {
    MetaData.TracerParticleDumpDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.TracerParticleDumpDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.HistoryDump.Name");
  if( strlen(dummy) > 0 ) {
    MetaData.HistoryDumpName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.HistoryDumpName, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.HistoryDump.Dir");
  if( strlen(dummy) > 0 ) {
    MetaData.HistoryDumpDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.HistoryDumpDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.LocalDir");
  if( strlen(dummy) > 0 ) {
    MetaData.LocalDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.LocalDir, dummy);
  }

  Param.GetScalar(dummy, "OutputControl.GlobalDir");
  if( strlen(dummy) > 0 ) {
    MetaData.GlobalDir = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.GlobalDir, dummy);
  }

  Param.GetScalar(LoadBalancing, "SimulationControl.Optimization.LoadBalancing");
  Param.GetScalar(ResetLoadBalancing, "SimulationControl.Optimization.ResetLoadBalancing");
  Param.GetScalar(LoadBalancingCycleSkip, "SimulationControl.Optimization.LoadBalancingCycleSkip");
  Param.GetScalar(LoadBalancingMinLevel, "SimulationControl.Optimization.LoadBalancingMinLevel");
  Param.GetScalar(LoadBalancingMaxLevel, "SimulationControl.Optimization.LoadBalancingMaxLevel");
  
  
  const int NumberOfTimeActions = Param.Size("SimulationControl.Timeaction.Actions");
  if (NumberOfTimeActions > MAX_TIME_ACTIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of TimeActions (%d > %d)!\n",NumberOfTimeActions,MAX_TIME_ACTIONS)
      }
  
  char TimeActionNames[MAX_LINE_LENGTH][MAX_TIME_ACTIONS];
  Param.GetArray(TimeActionNames, "SimulationControl.Timeaction.Actions");
  
  for (i = 0; i < NumberOfTimeActions; i++) {
    Param.GetScalar(TimeActionType[i], "SimulationControl.Timeaction.%s.Type", TimeActionNames[i]);
    Param.GetScalar(TimeActionRedshift[i], "SimulationControl.Timeaction.%s.Redshift", TimeActionNames[i]);
    Param.GetScalar(TimeActionTime[i], "SimulationControl.Timeaction.%s.Time", TimeActionNames[i]);
    Param.GetScalar(TimeActionParameter[i], "SimulationControl.Timeaction.%s.Parameter", TimeActionNames[i]);
  }
  
  Param.GetScalar(MetaData.StaticHierarchy, "SimulationControl.AMR.StaticHierarchy"); 
  
  Param.GetScalar(MetaData.TopGridRank, "SimulationControl.Domain.TopGridRank");
  Param.GetArray(MetaData.TopGridDims, "SimulationControl.Domain.TopGridDimensions");
  
  Param.GetScalar(MetaData.GravityBoundary, "Physics.TopGridGravityBoundary");
  
#ifdef TRANSFER
  Param.GetScalar(dummy, "Physics.RadiationField.RadHydroParamfile");
  if( strlen(dummy) > 0 ) {
    MetaData.RadHydroParameterFname = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RadHydroParameterFname, dummy);
  }
#endif

  Param.GetScalar(ImplicitProblem, "Physics.RadiationTransfer.ImplicitProblem");
  Param.GetScalar(RadiativeTransferFLD, "Physics.RadiationTransfer.RadiativeTransferFLD");
#ifdef EMISSIVITY
  Param.GetScalar(StarMakerEmissivityField, "Physics.OtherParticles.StarMaker.StarMakerEmissivityField");
  Param.GetScalar(uv_param, "Physics.RadiationField.uv_param");
#endif
  
  Param.GetScalar(MetaData.ParticleBoundaryType, "Physics.ParticleBoundaryType");
  Param.GetScalar(MetaData.NumberOfParticles, "Internal.NumberOfParticles");
  
  Param.GetScalar(MetaData.CourantSafetyNumber, "Physics.Hydro.CourantSafetyNumber");
  Param.GetScalar(MetaData.PPMFlatteningParameter, "Physics.Hydro.PPMFlatteningParameter"); 
  Param.GetScalar(MetaData.PPMDiffusionParameter, "Physics.Hydro.PPMDiffusionParameter"); 
  Param.GetScalar(MetaData.PPMSteepeningParameter, "Physics.Hydro.PPMSteepeningParameter"); 
  
  
  /* read global Parameters */
  Param.GetScalar(ProblemType, "Initialization.ProblemType");
  
#ifdef NEW_PROBLEM_TYPES
  if (sscanf(line, "ProblemTypeName = %s", dummy) == 1) {
    ProblemTypeName = dummy;
    ProblemType = -978;
    ret = 1;
  }
#endif
  
  Param.GetScalar(HydroMethod, "Physics.Hydro.HydroMethod");
  
  if (HydroMethod==MHD_RK) useMHD = 1;
  
  
  Param.GetScalar(huge_number, "SimulationControl.huge_number");
  Param.GetScalar(tiny_number, "SimulationControl.tiny_number");
  Param.GetScalar(Gamma, "Physics.Hydro.Gamma");
  Param.GetScalar(PressureFree, "Physics.Hydro.PressureFree");
  Param.GetScalar(RefineBy, "SimulationControl.AMR.RefineBy");
  Param.GetScalar(MaximumRefinementLevel, "SimulationControl.AMR.MaximumRefinementLevel");
  Param.GetScalar(MaximumGravityRefinementLevel, "SimulationControl.AMR.MaximumGravityRefinementLevel");
  Param.GetScalar(MaximumParticleRefinementLevel, "SimulationControl.AMR.MaximumParticleRefinementLevel");
  Param.GetArray(CellFlaggingMethod, "SimulationControl.AMR.CellFlaggingMethod");
  
  Param.GetScalar(FluxCorrection, "Physics.Hydro.FluxCorrection");
  Param.GetScalar(InterpolationMethod, "Physics.Hydro.InterpolationMethod");
  Param.GetScalar(ConservativeInterpolation, "Physics.Hydro.ConservativeInterpolation");
  Param.GetScalar(MinimumEfficiency, "SimulationControl.Optimization.MinimumEfficiency");
  Param.GetScalar(SubgridSizeAutoAdjust, "SimulationControl.Optimization.SubgridSizeAutoAdjust");
  Param.GetScalar(OptimalSubgridsPerProcessor, "SimulationControl.Optimization.OptimalSubgridsPerProcessor");
  Param.GetScalar(MinimumSubgridEdge, "SimulationControl.Optimization.MinimumSubgridEdge");
  Param.GetScalar(MaximumSubgridSize, "SimulationControl.Optimization.MaximumSubgridSize");
  Param.GetScalar(NumberOfBufferZones, "SimulationControl.AMR.NumberOfBufferZones");
  Param.GetScalar(FastSiblingLocatorEntireDomain, "SimulationControl.Optimization.FastSiblingLocatorEntireDomain");
  Param.GetScalar(MustRefineRegionMinRefinementLevel, "SimulationControl.AMR.MustRefineRegionMinRefinementLevel");
  Param.GetScalar(MetallicityRefinementMinLevel, "SimulationControl.AMR.MetallicityRefinementMinLevel");  
  Param.GetScalar(MetallicityRefinementMinMetallicity, "SimulationControl.AMR.MetallicityRefinementMinMetallicity"); 
  Param.GetScalar(MetallicityRefinementMinDensity, "SimulationControl.AMR.MetallicityRefinementMinDensity"); 
  Param.GetArray(DomainLeftEdge, "SimulationControl.Domain.DomainLeftEdge");
  Param.GetArray(DomainRightEdge, "SimulationControl.Domain.DomainRightEdge");
  Param.GetArray(GridVelocity, "Physics.GridVelocity");
  
  Param.GetScalar(RefineRegionAutoAdjust, "SimulationControl.Domain.RefineRegionAutoAdjust"); 
  Param.GetArray(RefineRegionLeftEdge, "SimulationControl.Domain.RefineRegionLeftEdge");
  Param.GetArray(RefineRegionRightEdge, "SimulationControl.Domain.RefineRegionRightEdge");
  Param.GetArray(MustRefineRegionLeftEdge, "SimulationControl.Domain.MustRefineRegionLeftEdge");
  Param.GetArray(MustRefineRegionRightEdge, "SimulationControl.Domain.MustRefineRegionRightEdge");

  
  /* Read evolving RefineRegion */
  
  Param.GetScalar(RefineRegionTimeType, "SimulationControl.AMR.RefineRegionTimeType");
    
  Param.GetScalar(dummy, "SimulationControl.AMR.RefineRegionFile");
  if( strlen(dummy) > 0 ) {
    RefineRegionFile = new char[MAX_LINE_LENGTH];
    strcpy(RefineRegionFile, dummy);
  }

  int NumberOfFields = Param.Size("Internal.Fields");
  if (NumberOfFields > MAX_NUMBER_OF_BARYON_FIELDS) {
    ENZO_VFAIL("You've exceeded the maximum number of BaryonFields (%d > %d)!\n",NumberOfFields,MAX_NUMBER_OF_BARYON_FIELDS)
      }
  
  char FieldNames[MAX_LINE_LENGTH][MAX_NUMBER_OF_BARYON_FIELDS];
  Param.GetArray(FieldNames, "Internal.Fields");
  
  // for (i = 0; i < NumberOfFields; i++) {
  //   printf("Found field : %s\n",FieldNames[i]);
  // }
	
  for (i = 0; i < NumberOfFields; i++) {
    Param.GetScalar(dummy, "Internal.%s.Name", FieldNames[i]);
    if( strlen(dummy) > 0 ) {
      DataLabel[i] = new char[MAX_LINE_LENGTH];
      strcpy(DataLabel[i], dummy);
    }

    Param.GetScalar(dummy, "Internal.%s.cgsConversionFactor", FieldNames[i]);
    if( strlen(dummy) > 0 ) {
      DataUnits[i] = new char[MAX_LINE_LENGTH];
      strcpy(DataUnits[i], dummy);
    }
  }
  
  Param.GetScalar(UniformGravity, "Physics.Gravity.UniformGravity");
  Param.GetScalar(UniformGravityDirection, "Physics.Gravity.UniformGravityDirection");
  Param.GetScalar(UniformGravityConstant, "Physics.Gravity.UniformGravityConstant");
  
  Param.GetScalar(PointSourceGravity, "Physics.Gravity.PointSourceGravity");
  Param.GetArray(PointSourceGravityPosition, "Physics.Gravity.PointSourceGravityPosition");
  Param.GetScalar(PointSourceGravityConstant, "Physics.Gravity.PointSourceGravityConstant");
  Param.GetScalar(PointSourceGravityCoreRadius, "Physics.Gravity.PointSourceGravityCoreRadius");
  
  Param.GetScalar(ExternalGravity, "Physics.Gravity.ExternalGravity");
  
  Param.GetScalar(SelfGravity, "Physics.Gravity.SelfGravity");
  Param.GetScalar(SelfGravityGasOff, "Physics.Gravity.SelfGravityGasOff");
  Param.GetScalar(AccretionKernel, "Physics.AccretionKernel");
  Param.GetScalar(GravitationalConstant, "Physics.Gravity.GravitationalConstant");
  
  Param.GetScalar(S2ParticleSize, "Physics.OtherParticles.S2ParticleSize");
  Param.GetScalar(GravityResolution, "Physics.Gravity.GravityResolution");
  Param.GetScalar(ComputePotential, "Physics.Gravity.ComputePotential");
  Param.GetScalar(PotentialIterations, "Physics.Gravity.PotentialIterations");
  Param.GetScalar(WritePotential, "OutputControl.SupplementalFields.WritePotential"); 
  Param.GetScalar(BaryonSelfGravityApproximation, "Physics.Gravity.BaryonSelfGravityApproximation"); 
  
  Param.GetScalar(GreensFunctionMaxNumber, "Physics.Gravity.GreensFunctionMaxNumber");
  Param.GetScalar(GreensFunctionMaxSize, "Physics.Gravity.GreensFunctionMaxSize");
  
  Param.GetScalar(DualEnergyFormalism, "Physics.Hydro.DualEnergyFormalism");
  Param.GetScalar(DualEnergyFormalismEta1, "Physics.Hydro.DualEnergyFormalismEta1");
  Param.GetScalar(DualEnergyFormalismEta2, "Physics.Hydro.DualEnergyFormalismEta2");
  Param.GetScalar(ParticleCourantSafetyNumber, "Physics.Hydro.ParticleCourantSafetyNumber");
  Param.GetScalar(RootGridCourantSafetyNumber, "Physics.Hydro.RootGridCourantSafetyNumber");
  Param.GetScalar(RandomForcing, "Physics.Miscellaneous.RandomForcing");
  Param.GetScalar(RandomForcingEdot, "Physics.Miscellaneous.RandomForcingEdot");
  Param.GetScalar(RandomForcingMachNumber, "Physics.Miscellaneous.RandomForcingMachNumber");
  Param.GetScalar(RadiativeCooling, "Physics.AtomicPhysics.RadiativeCooling");
  Param.GetScalar(RadiativeCoolingModel, "Physics.AtomicPhysics.RadiativeCoolingModel");
  Param.GetScalar(GadgetEquilibriumCooling, "Physics.AtomicPhysics.GadgetEquilibriumCooling");
  Param.GetScalar(MultiSpecies, "Physics.AtomicPhysics.MultiSpecies");
  Param.GetScalar(PrimordialChemistrySolver, "Physics.AtomicPhysics.PrimordialChemistrySolver");
  Param.GetScalar(CIECooling, "Physics.AtomicPhysics.CIECooling");
  Param.GetScalar(H2OpticalDepthApproximation, "Physics.AtomicPhysics.H2OpticalDepthApproximation");
  Param.GetScalar(ThreeBodyRate, "Physics.AtomicPhysics.ThreeBodyRate");
  

  Param.GetScalar(dummy, "Physics.AtomicPhysics.CloudyCooling.CloudyCoolingGridFile");
  if( strlen(dummy) > 0 ) {
    CloudyCoolingData.CloudyCoolingGridFile = new char[MAX_LINE_LENGTH];
    strcpy(CloudyCoolingData.CloudyCoolingGridFile, dummy);
  }

  Param.GetScalar(CloudyCoolingData.IncludeCloudyHeating, "Physics.AtomicPhysics.CloudyCooling.IncludeCloudyHeating");
  Param.GetScalar(CloudyCoolingData.IncludeCloudyMMW, "Physics.AtomicPhysics.CloudyCooling.IncludeCloudyMMW");
  Param.GetScalar(CloudyCoolingData.CMBTemperatureFloor, "Physics.AtomicPhysics.CloudyCooling.CMBTemperatureFloor");
  Param.GetScalar(CloudyCoolingData.CloudyMetallicityNormalization, "Physics.AtomicPhysics.CloudyCooling.CloudyMetallicityNormalization");
  Param.GetScalar(CloudyCoolingData.CloudyElectronFractionFactor, "Physics.AtomicPhysics.CloudyCooling.CloudyElectronFractionFactor");
  Param.GetScalar(MetalCooling, "Physics.AtomicPhysics.MetalCooling");

  Param.GetScalar(dummy, "Physics.AtomicPhysics.MetalCoolingTable");
  if( strlen(dummy) > 0 ) {
    MetalCoolingTable = new char[MAX_LINE_LENGTH];
    strcpy(MetalCoolingTable, dummy);
  }

  Param.GetScalar(CRModel, "Physics.Miscellaneous.CRModel");
  Param.GetScalar(ShockMethod, "Physics.Miscellaneous.ShockMethod");
  Param.GetScalar(ShockTemperatureFloor, "Physics.Miscellaneous.ShockTemperatureFloor");
  Param.GetScalar(StorePreShockFields, "Physics.Miscellaneous.StorePreShockFields");

  Param.GetScalar(RadiationFieldType, "Physics.RadiationField.RadiationFieldType");
  Param.GetScalar(TabulatedLWBackground, "Physics.RadiationField.TabulatedLWBackground");
  Param.GetScalar(AdjustUVBackground, "Physics.RadiationField.AdjustUVBackground");
  Param.GetScalar(SetUVBAmplitude, "Physics.RadiationField.SetUVBAmplitude");
  Param.GetScalar(SetHeIIHeatingScale, "Physics.RadiationField.SetHeIIHeatingScale");
  Param.GetScalar(RadiationFieldLevelRecompute, "Physics.RadiationField.RadiationFieldLevelRecompute");
  Param.GetScalar(RadiationData.RadiationShield, "Physics.RadiationField.RadiationShield");
  Param.GetScalar(CoolData.f3, "Physics.RadiationField.RadiationSpectrumNormalization");
  Param.GetScalar(CoolData.alpha0, "Physics.RadiationField.RadiationSpectrumSlope");
  Param.GetScalar(CoolData.f0to3, "Physics.AtomicPhysics.CoolDataf0to3");
  Param.GetScalar(CoolData.RadiationRedshiftOn, "Physics.RadiationField.RadiationRedshiftOn");
  Param.GetScalar(CoolData.RadiationRedshiftOff, "Physics.RadiationField.RadiationRedshiftOff");
  Param.GetScalar(CoolData.RadiationRedshiftFullOn, "Physics.RadiationField.RadiationRedshiftFullOn");
  Param.GetScalar(CoolData.RadiationRedshiftDropOff, "Physics.RadiationField.RadiationRedshiftDropOff");
  Param.GetScalar(CoolData.HydrogenFractionByMass, "Physics.AtomicPhysics.HydrogenFractionByMass");
  Param.GetScalar(CoolData.DeuteriumToHydrogenRatio, "Physics.AtomicPhysics.DeuteriumToHydrogenRatio");
  Param.GetScalar(CoolData.NumberOfTemperatureBins, "Physics.AtomicPhysics.NumberOfTemperatureBins");
  Param.GetScalar(CoolData.ih2co, "Physics.AtomicPhysics.CoolDataIh2co");
  Param.GetScalar(CoolData.ipiht, "Physics.AtomicPhysics.CoolDataIpiht");
  Param.GetScalar(CoolData.TemperatureStart, "Physics.AtomicPhysics.TemperatureStart");
  Param.GetScalar(CoolData.TemperatureEnd, "Physics.AtomicPhysics.TemperatureEnd");
  Param.GetScalar(CoolData.comp_xray, "Physics.AtomicPhysics.CoolDataCompXray");
  Param.GetScalar(CoolData.temp_xray, "Physics.AtomicPhysics.CoolDataTempXray");
  Param.GetScalar(RateData.CaseBRecombination, "Physics.AtomicPhysics.RateDataCaseBRecombination");
  Param.GetScalar(PhotoelectricHeating, "Physics.AtomicPhysics.PhotoelectricHeating");
  
  Param.GetScalar(OutputCoolingTime, "OutputControl.SupplementalFields.OutputCoolingTime"); 
  Param.GetScalar(OutputTemperature, "OutputControl.SupplementalFields.OutputTemperature"); 
  Param.GetScalar(OutputSmoothedDarkMatter, "OutputControl.SupplementalFields.OutputSmoothedDarkMatter"); 
  Param.GetScalar(SmoothedDarkMatterNeighbors, "OutputControl.SupplementalFields.SmoothedDarkMatterNeighbors"); 
  Param.GetScalar(OutputGriddedStarParticle, "OutputControl.SupplementalFields.OutputGriddedStarParticle"); 
  
  Param.GetScalar(ZEUSQuadraticArtificialViscosity, "Physics.Hydro.ZEUSQuadraticArtificialViscosity");
  Param.GetScalar(ZEUSLinearArtificialViscosity, "Physics.Hydro.ZEUSLinearArtificialViscosity");

  Param.GetScalar(UseMinimumPressureSupport, "Physics.Hydro.UseMinimumPressureSupport"); 
  Param.GetScalar(MinimumPressureSupportParameter, "Physics.Hydro.MinimumPressureSupportParameter");
  Param.GetScalar(RefineByJeansLengthSafetyFactor, "SimulationControl.AMR.RefineByJeansLengthSafetyFactor");
  Param.GetScalar(RefineByJeansLengthUnits, "SimulationControl.AMR.RefineByJeansLengthUnits"); 
  Param.GetScalar(JeansRefinementColdTemperature, "SimulationControl.AMR.JeansRefinementColdTemperature");
  Param.GetScalar(RefineByResistiveLengthSafetyFactor, "SimulationControl.AMR.RefineByResistiveLengthSafetyFactor");
  Param.GetScalar(MustRefineParticlesRefineToLevel, "SimulationControl.AMR.MustRefineParticlesRefineToLevel");
  Param.GetScalar(MustRefineParticlesRefineToLevelAutoAdjust, "SimulationControl.AMR.MustRefineParticlesRefineToLevelAutoAdjust"); 
  Param.GetScalar(MustRefineParticlesMinimumMass, "SimulationControl.AMR.MustRefineParticlesMinimumMass");
  Param.GetScalar(ParticleTypeInFile, "OutputControl.ParticleTypeInFile"); 

  int NumberOfStaticRefineRegions = Param.Size("SimulationControl.AMR.StaticRefineRegion.Regions");
  if (NumberOfStaticRefineRegions > MAX_STATIC_REGIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of StaticRefineRegions (%d > %d)!\n", NumberOfStaticRefineRegions, MAX_STATIC_REGIONS)
      }
  
  char StaticRefineRegionNames[MAX_LINE_LENGTH][MAX_STATIC_REGIONS];//NumberOfStaticRefineRegions];
  Param.GetArray(StaticRefineRegionNames, "SimulationControl.AMR.StaticRefineRegion.Regions");
  
  for (i = 0; i < NumberOfStaticRefineRegions; i++) {
    Param.GetScalar(StaticRefineRegionLevel[i], "SimulationControl.AMR.StaticRefineRegion.%s.Level",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionLeftEdge[i],"SimulationControl.AMR.StaticRefineRegion.%s.LeftEdge",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionRightEdge[i],"SimulationControl.AMR.StaticRefineRegion.%s.RightEdge",StaticRefineRegionNames[i]);
  }
  
  Param.GetScalar(ParallelRootGridIO, "SimulationControl.Optimization.ParallelRootGridIO");
  Param.GetScalar(ParallelParticleIO, "SimulationControl.Optimization.ParallelParticleIO");
  
  Param.GetScalar(Unigrid, "SimulationControl.Optimization.Unigrid");
  Param.GetScalar(UnigridTranspose, "SimulationControl.Optimization.UnigridTranspose");
  Param.GetScalar(NumberOfRootGridTilesPerDimensionPerProcessor, "SimulationControl.Optimization.NumberOfRootGridTilesPerDimensionPerProcessor");
  
  Param.GetScalar(PartitionNestedGrids, "Initialization.PartitionNestedGrids");
  Param.GetScalar(StaticPartitionNestedGrids, "Initialization.StaticPartitionNestedGrids");
  
  Param.GetScalar(ExtractFieldsOnly, "OutputControl.ExtractFieldsOnly");
  
  Param.GetScalar(debug1, "Internal.Debug1");
  
  Param.GetScalar(debug2, "Internal.Debug2");
  
  Param.GetScalar(MemoryLimit, "SimulationControl.Optimization.MemoryLimit");
  
  
#ifdef OOC_BOUNDARY
  
  Param.GetScalar(ExternalBoundaryIO, "OutputControl.SupplementalFields.ExternalBoundaryIO");
  
  Param.GetScalar(ExternalBoundaryTypeIO, "OutputControl.SupplementalFields.ExternalBoundaryTypeIO");
  
  Param.GetScalar(ExternalBoundaryValueIO, "OutputControl.SupplementalFields.ExternalBoundaryValueIO");
  
  Param.GetScalar(SimpleConstantBoundary, "SimulationControl.Domain.SimpleConstantBoundary");
  
#endif
  
  Param.GetArray(SlopeFlaggingFields, "SimulationControl.AMR.SlopeFlaggingFields");
  
  Param.GetArray(MinimumSlopeForRefinement, "SimulationControl.AMR.MinimumSlopeForRefinement");
  
  Param.GetArray(MinimumOverDensityForRefinement, "SimulationControl.AMR.MinimumOverDensityForRefinement");
  
  Param.GetArray(MinimumMassForRefinement, "SimulationControl.AMR.MinimumMassForRefinement");
  
  Param.GetArray(MinimumMassForRefinementLevelExponent, "SimulationControl.AMR.MinimumMassForRefinementLevelExponent");
  
  
  Param.GetScalar(MinimumPressureJumpForRefinement, "SimulationControl.AMR.MinimumPressureJumpForRefinement");
  Param.GetScalar(MinimumShearForRefinement, "SimulationControl.AMR.MinimumShearForRefinement");
  Param.GetScalar(MinimumEnergyRatioForRefinement, "SimulationControl.AMR.MinimumEnergyRatioForRefinement");
  Param.GetScalar(ShockwaveRefinementMinMach, "SimulationControl.AMR.ShockwaveRefinementMinMach");
  Param.GetScalar(ShockwaveRefinementMinVelocity, "SimulationControl.AMR.ShockwaveRefinementMinVelocity");
  Param.GetScalar(ShockwaveRefinementMaxLevel, "SimulationControl.AMR.ShockwaveRefinementMaxLevel");
  Param.GetScalar(ComovingCoordinates, "Physics.Cosmology.ComovingCoordinates");
  Param.GetScalar(StarParticleCreation, "Physics.OtherParticles.StarParticleCreation");
  Param.GetScalar(BigStarFormation, "Physics.OtherParticles.BigStar.Formation");
  Param.GetScalar(BigStarFormationDone, "Physics.OtherParticles.BigStar.FormationDone");
  Param.GetScalar(BigStarSeparation, "Physics.OtherParticles.BigStar.Separation");
  Param.GetScalar(SimpleQ, "Physics.RadiationTransfer.SimpleQ");
  Param.GetScalar(SimpleRampTime, "Physics.RadiationTransfer.SimpleRampTime");
  Param.GetScalar(StarParticleFeedback, "Physics.OtherParticles.StarParticleFeedback");
  Param.GetScalar(NumberOfParticleAttributes, "Physics.OtherParticles.NumberOfParticleAttributes");
  Param.GetScalar(AddParticleAttributes, "Physics.OtherParticles.AddParticleAttributes");
  
  
  /* read data which defines the boundary conditions */
  Param.GetArray(MetaData.LeftFaceBoundaryCondition, "SimulationControl.Domain.LeftFaceBoundaryCondition");
  Param.GetArray(MetaData.RightFaceBoundaryCondition, "SimulationControl.Domain.RightFaceBoundaryCondition");

  Param.GetScalar(dummy, "Internal.BoundaryConditionName");
  if( strlen(dummy) > 0 ) {
    MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.BoundaryConditionName, dummy);
  }
	
  Param.GetScalar(dummy, "Internal.Provenance.MetaDataIdentifier");
  if( strlen(dummy) > 0 ) {
    MetaData.MetaDataIdentifier = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.MetaDataIdentifier, dummy);
  }

  Param.GetScalar(dummy, "Internal.Provenance.SimulationUUID");
  if( strlen(dummy) > 0 ) {
    MetaData.SimulationUUID = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.SimulationUUID, dummy);
  }

  // Param.GetScalar(dummy, "Internal.Provenance.DatasetUUID");
  // if( strlen(dummy) > 0 ) {
  //   MetaData.DatasetUUID = new char[MAX_LINE_LENGTH];
  //   strcpy(MetaData.DatasetUUID, dummy);
  // }

  Param.GetScalar(dummy, "Internal.Provenance.RestartDatasetUUID");
  if( strlen(dummy) > 0 ) {
    MetaData.RestartDatasetUUID = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RestartDatasetUUID, dummy);
  }

  Param.GetScalar(dummy, "Internal.Provenance.InitialConditionsUUID");
  if( strlen(dummy) > 0 ) {
    MetaData.InitialConditionsUUID = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.InitialConditionsUUID, dummy);
  }

  /* Check version number. */
  
  Param.GetScalar(TempFloat, "Internal.Provenance.VersionNumber");
  if (fabs(TempFloat - VERSION) >= 1.0e-3 &&
      MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Warning: Incorrect version number.\n");
  
  
  /* Read star particle parameters. */
  
  Param.GetScalar(StarMakerOverDensityThreshold, "Physics.OtherParticles.StarMaker.OverDensityThreshold");
  Param.GetScalar(StarMakerSHDensityThreshold, "Physics.OtherParticles.StarMaker.SHDensityThreshold");
  Param.GetScalar(StarMakerMassEfficiency, "Physics.OtherParticles.StarMaker.MassEfficiency");
  Param.GetScalar(StarMakerMinimumMass, "Physics.OtherParticles.StarMaker.MinimumMass");
  Param.GetScalar(StarMakerMinimumDynamicalTime, "Physics.OtherParticles.StarMaker.MinimumDynamicalTime");
  Param.GetScalar(StarMassEjectionFraction, "Physics.OtherParticles.StarMassEjectionFraction");
  Param.GetScalar(StarMetalYield, "Physics.OtherParticles.StarMetalYield");
  Param.GetScalar(StarEnergyToThermalFeedback, "Physics.OtherParticles.StarEnergyToThermalFeedback");
  Param.GetScalar(StarEnergyToStellarUV, "Physics.OtherParticles.StarEnergyToStellarUV");
  Param.GetScalar(StarEnergyToQuasarUV, "Physics.OtherParticles.StarEnergyToQuasarUV");
  
  Param.GetScalar(StarFeedbackDistRadius, "Physics.OtherParticles.StarFeedbackDistRadius");
  Param.GetScalar(StarFeedbackDistCellStep, "Physics.OtherParticles.StarFeedbackDistCellStep");
  Param.GetScalar(StarClusterUseMetalField, "Physics.OtherParticles.StarClusterParticle.UseMetalField");
  Param.GetScalar(StarClusterMinDynamicalTime, "Physics.OtherParticles.StarClusterParticle.MinDynamicalTime");
  Param.GetScalar(StarClusterHeliumIonization, "Physics.OtherParticles.StarClusterParticle.HeliumIonization");
  Param.GetScalar(StarClusterUnresolvedModel, "Physics.OtherParticles.StarClusterParticle.UnresolvedModel");
  Param.GetScalar(StarClusterIonizingLuminosity, "Physics.OtherParticles.StarClusterParticle.IonizingLuminosity");
  Param.GetScalar(StarClusterSNEnergy, "Physics.OtherParticles.StarClusterParticle.SNEnergy");
  Param.GetScalar(StarClusterSNRadius, "Physics.OtherParticles.StarClusterParticle.SNRadius");
  Param.GetScalar(StarClusterFormEfficiency, "Physics.OtherParticles.StarClusterParticle.FormEfficiency");
  Param.GetScalar(StarClusterMinimumMass, "Physics.OtherParticles.StarClusterParticle.MinimumMass");
  Param.GetScalar(StarClusterCombineRadius, "Physics.OtherParticles.StarClusterParticle.CombineRadius");
  
  Param.GetArray(StarClusterRegionLeftEdge, "Physics.OtherParticles.StarClusterParticle.RegionLeftEdge");
  Param.GetArray(StarClusterRegionRightEdge, "Physics.OtherParticles.StarClusterParticle.RegionRightEdge");
  
  Param.GetScalar(PopIIIStarMass, "Physics.OtherParticles.PopIIIStarParticle.StarMass");                           
  Param.GetScalar(PopIIIInitialMassFunctionSeed, "Physics.OtherParticles.PopIIIStarParticle.InitialMassFunctionSeed");                 
  Param.GetScalar(PopIIIInitialMassFunctionCalls, "Physics.OtherParticles.PopIIIStarParticle.InitialMassFunctionCalls"); 
  Param.GetScalar(PopIIILowerMassCutoff, "Physics.OtherParticles.PopIIIStarParticle.LowerMassCutoff");
  Param.GetScalar(PopIIIUpperMassCutoff, "Physics.OtherParticles.PopIIIStarParticle.UpperMassCutoff");
  
  //   ret += sscanf(line,"PopIIIMassRange %"FSYM" %"FSYM, &PopIIILowerMassCutoff, &UpperMassCutoff);
  
  Param.GetScalar(PopIIIInitialMassFunctionSlope, "Physics.OtherParticles.PopIIIStarParticle.InitialMassFunctionSlope");
  Param.GetScalar(PopIIIBlackHoles, "Physics.OtherParticles.PopIIIStarParticle.BlackHoles");
  Param.GetScalar(PopIIIBHLuminosityEfficiency, "Physics.OtherParticles.PopIIIStarParticle.BHLuminosityEfficiency");
  Param.GetScalar(PopIIIOverDensityThreshold, "Physics.OtherParticles.PopIIIStarParticle.OverDensityThreshold");
  Param.GetScalar(PopIIIH2CriticalFraction, "Physics.OtherParticles.PopIIIStarParticle.H2CriticalFraction");
  Param.GetScalar(PopIIIMetalCriticalFraction, "Physics.OtherParticles.PopIIIStarParticle.MetalCriticalFraction");
  Param.GetScalar(PopIIISupernovaRadius, "Physics.OtherParticles.PopIIIStarParticle.SupernovaRadius");
  Param.GetScalar(PopIIISupernovaUseColour, "Physics.OtherParticles.PopIIIStarParticle.SupernovaUseColour");
  Param.GetScalar(PopIIISupernovaMustRefine, "Physics.OtherParticles.PopIIIStarParticle.SupernovaMustRefine");
  Param.GetScalar(PopIIISupernovaMustRefineResolution, "Physics.OtherParticles.PopIIIStarParticle.SupernovaMustRefineResolution");
  Param.GetScalar(PopIIIHeliumIonization, "Physics.OtherParticles.PopIIIStarParticle.HeliumIonization");
  Param.GetScalar(PopIIIColorDensityThreshold, "Physics.OtherParticles.PopIIIStarParticle.ColorDensityThreshold");
  Param.GetScalar(PopIIIColorMass, "Physics.OtherParticles.PopIIIStarParticle.ColorMass");
  
  Param.GetScalar(MBHAccretion, "Physics.OtherParticles.MBHParticle.Accretion");
  Param.GetScalar(MBHAccretionRadius, "Physics.OtherParticles.MBHParticle.AccretionRadius");
  Param.GetScalar(MBHAccretingMassRatio, "Physics.OtherParticles.MBHParticle.AccretingMassRatio");
  Param.GetScalar(MBHAccretionFixedTemperature, "Physics.OtherParticles.MBHParticle.AccretionFixedTemperature");
  Param.GetScalar(MBHAccretionFixedRate, "Physics.OtherParticles.MBHParticle.AccretionFixedRate");
  Param.GetScalar(MBHTurnOffStarFormation, "Physics.OtherParticles.MBHParticle.TurnOffStarFormation");
  Param.GetScalar(MBHCombineRadius, "Physics.OtherParticles.MBHParticle.CombineRadius");
  Param.GetScalar(MBHMinDynamicalTime, "Physics.OtherParticles.MBHParticle.MinDynamicalTime");
  Param.GetScalar(MBHMinimumMass, "Physics.OtherParticles.MBHParticle.MinimumMass");
  
  Param.GetScalar(MBHFeedback, "Physics.OtherParticles.MBHParticle.Feedback");
  Param.GetScalar(MBHFeedbackRadiativeEfficiency, "Physics.OtherParticles.MBHParticle.FeedbackRadiativeEfficiency");
  Param.GetScalar(MBHFeedbackEnergyCoupling, "Physics.OtherParticles.MBHParticle.FeedbackEnergyCoupling");
  Param.GetScalar(MBHFeedbackMassEjectionFraction, "Physics.OtherParticles.MBHParticle.FeedbackMassEjectionFraction");
  Param.GetScalar(MBHFeedbackMetalYield, "Physics.OtherParticles.MBHParticle.FeedbackMetalYield");
  Param.GetScalar(MBHFeedbackThermalRadius, "Physics.OtherParticles.MBHParticle.FeedbackThermalRadius");
  Param.GetScalar(MBHFeedbackJetsThresholdMass, "Physics.OtherParticles.MBHParticle.FeedbackJetsThresholdMass");
  Param.GetScalar(MBHParticleIO, "Initialization.MBHParticleIO");
  
  Param.GetScalar(dummy, "Initialization.MBHParticleIOFilename");
  if( strlen(dummy) > 0 ) {
    MBHParticleIOFilename = new char[MAX_LINE_LENGTH];
    strcpy(MBHParticleIOFilename, dummy);
  }

  Param.GetScalar(dummy, "Initialization.MBHInsertLocationFilename");
  if( strlen(dummy) > 0 ) {
    MBHInsertLocationFilename = new char[MAX_LINE_LENGTH];
    strcpy(MBHInsertLocationFilename, dummy);
  }

  
  
  /* Read Movie Dump parameters */
  
  Param.GetScalar(MovieSkipTimestep, "OutputControl.MovieDump.SkipTimestep");
  Param.GetScalar(Movie3DVolumes, "OutputControl.MovieDump.Volumes3D");
  Param.GetScalar(MovieVertexCentered, "OutputControl.MovieDump.VertexCentered");
  Param.GetScalar(NewMovieParticleOn, "OutputControl.MovieDump.ParticleOn");
  Param.GetArray(MovieDataField, "OutputControl.MovieDump.DataField");
  Param.GetScalar(NewMovieDumpNumber, "OutputControl.MovieDump.DumpNumber");

  Param.GetScalar(dummy, "OutputControl.MovieDump.Name");
  if( strlen(dummy) > 0 ) {
    NewMovieName = new char[MAX_LINE_LENGTH];
    strcpy(NewMovieName, dummy);
  }

  Param.GetScalar(MetaData.MovieTimestepCounter, "Internal.OutputLabeling.MovieTimestepCounter");
  
  Param.GetScalar(MultiMetals, "Physics.AtomicPhysics.MultiMetals");
  Param.GetScalar(IsotropicConduction, "Physics.Conduction.IsotropicConduction");
  Param.GetScalar(IsotropicConductionSpitzerFraction, "Physics.Conduction.IsotropicConductionSpitzerFraction");
  Param.GetScalar(AnisotropicConduction, "Physics.Conduction.AnisotropicConduction");
  Param.GetScalar(AnisotropicConductionSpitzerFraction, "Physics.Conduction.AnisotropicConductionSpitzerFraction");
  Param.GetScalar(ConductionCourantSafetyNumber, "Physics.Conduction.ConductionCourantSafetyNumber");
  
  Param.GetScalar(RadiativeTransfer, "Physics.RadiationTransfer.RadiativeTransfer");
  Param.GetScalar(RadiationXRaySecondaryIon, "Physics.RadiationTransfer.RadiationXRaySecondaryIon");
  Param.GetScalar(RadiationXRayComptonHeating, "Physics.RadiationTransfer.RadiationXRayComptonHeating");
  
  
  /* Shearing Box Boundary parameters */
  
  Param.GetScalar(AngularVelocity, "Physics.AngularVelocity");
  Param.GetScalar(VelocityGradient, "Physics.VelocityGradient");
  Param.GetScalar(ShearingVelocityDirection, "Physics.ShearingVelocityDirection");
  Param.GetScalar(ShearingBoxProblemType, "Physics.ShearingBoxProblemType");
  
  /* Embedded Python */
  Param.GetScalar(PythonTopGridSkip, "Analysis.Python.PythonTopGridSkip");
  Param.GetScalar(PythonSubcycleSkip, "Analysis.Python.PythonSubcycleSkip");

#ifdef USE_PYTHON
  Param.GetScalar(NumberOfPythonCalls, "Analysis.Python.NumberOfPythonCalls");
  Param.GetScalar(NumberOfPythonTopGridCalls, "Analysis.Python.NumberOfPythonTopGridCalls");
  Param.GetScalar(NumberOfPythonSubcycleCalls, "Analysis.Python.NumberOfPythonSubcycleCalls");
#endif

    /* Inline halo finder */

  Param.GetScalar(InlineHaloFinder, "Analysis.HaloFinder.InlineHaloFinder");
  Param.GetScalar(HaloFinderSubfind, "Analysis.HaloFinder.HaloFinderSubfind");
  Param.GetScalar(HaloFinderOutputParticleList, "Analysis.HaloFinder.HaloFinderOutputParticleList");
  Param.GetScalar(HaloFinderRunAfterOutput, "Analysis.HaloFinder.HaloFinderRunAfterOutput");
  Param.GetScalar(HaloFinderLinkingLength, "Analysis.HaloFinder.HaloFinderLinkingLength");
  Param.GetScalar(HaloFinderMinimumSize, "Analysis.HaloFinder.HaloFinderMinimumSize");
  Param.GetScalar(HaloFinderCycleSkip, "Analysis.HaloFinder.HaloFinderCycleSkip");
  Param.GetScalar(HaloFinderTimestep, "Analysis.HaloFinder.HaloFinderTimestep");
  Param.GetScalar(HaloFinderLastTime, "Analysis.HaloFinder.HaloFinderLastTime");
  
  /* This Block for Stanford Hydro */
  
  Param.GetScalar(UseHydro, "Physics.Hydro.UseHydro");
  
  
  /* Sink particles (for present day star formation) & winds */
  Param.GetScalar(SinkMergeDistance, "Physics.OtherParticles.SinkParticle.MergeDistance");
  Param.GetScalar(SinkMergeMass, "Physics.OtherParticles.SinkParticle.MergeMass");
  Param.GetScalar(StellarWindFeedback, "Physics.OtherParticles.StellarWindFeedback");
  Param.GetScalar(StellarWindTurnOnMass, "Physics.OtherParticles.StellarWindTurnOnMass");
  Param.GetScalar(MSStellarWindTurnOnMass, "Physics.OtherParticles.MSStellarWindTurnOnMass");
  
  Param.GetScalar(VelAnyl, "OutputControl.SupplementalFields.VelAnyl");
  Param.GetScalar(BAnyl, "OutputControl.SupplementalFields.BAnyl");
  
  
  
  /* Read MHD Paramters */
  Param.GetScalar(UseDivergenceCleaning, "Physics.MHD.UseDivergenceCleaning");
  Param.GetScalar(DivergenceCleaningBoundaryBuffer, "Physics.MHD.DivergenceCleaningBoundaryBuffer");
  Param.GetScalar(DivergenceCleaningThreshold, "Physics.MHD.DivergenceCleaningThreshold");
  Param.GetScalar(PoissonApproximationThreshold, "Physics.PoissonApproximationThreshold");
  Param.GetScalar(PoissonBoundaryType, "Physics.PoissonBoundaryType");
  
  
  Param.GetScalar(AngularVelocity, "Physics.AngularVelocity");
  Param.GetScalar(VelocityGradient, "Physics.VelocityGradient");
  Param.GetScalar(UseDrivingField, "Physics.Miscellaneous.UseDrivingField");
  Param.GetScalar(DrivingEfficiency, "Physics.Miscellaneous.DrivingEfficiency");
  
  Param.GetScalar(StringKick, "Initialization.StringKick");
  Param.GetScalar(StringKickDimension, "Initialization.StringKickDimension");
  Param.GetScalar(UsePhysicalUnit, "Physics.UsePhysicalUnit");
  Param.GetScalar(Theta_Limiter, "Physics.Hydro.Theta_Limiter");
  Param.GetScalar(RKOrder, "Physics.Hydro.RKOrder");
  Param.GetScalar(UseFloor, "Physics.UseFloor");
  Param.GetScalar(UseViscosity, "Physics.Hydro.UseViscosity");
  Param.GetScalar(ViscosityCoefficient, "Physics.Hydro.ViscosityCoefficient");
  
  Param.GetScalar(UseAmbipolarDiffusion, "Physics.Hydro.UseAmbipolarDiffusion");
  Param.GetScalar(UseResistivity, "Physics.Hydro.UseResistivity");
  Param.GetScalar(SmallRho, "Physics.SmallRho");
  Param.GetScalar(SmallP, "Physics.SmallP");
  Param.GetScalar(SmallT, "Physics.SmallT");
  Param.GetScalar(MaximumAlvenSpeed, "Physics.MHD.MaximumAlvenSpeed");
  Param.GetScalar(Coordinate, "Physics.Coordinate");
  Param.GetScalar(RiemannSolver, "Physics.Hydro.RiemannSolver");
  Param.GetScalar(RiemannSolverFallback, "Physics.Hydro.RiemannSolverFallback");
  Param.GetScalar(ConservativeReconstruction, "Physics.Hydro.ConservativeReconstruction");
  Param.GetScalar(PositiveReconstruction, "Physics.Hydro.PositiveReconstruction");
  Param.GetScalar(ReconstructionMethod, "Physics.Hydro.ReconstructionMethod");
  
  Param.GetScalar(EOSType, "Physics.Hydro.EOSType");
  Param.GetScalar(EOSSoundSpeed, "Physics.Hydro.EOSSoundSpeed");
  Param.GetScalar(EOSCriticalDensity, "Physics.Hydro.EOSCriticalDensity");
  Param.GetScalar(EOSGamma, "Physics.Hydro.EOSGamma");
  Param.GetScalar(UseConstantAcceleration, "Physics.Gravity.UseConstantAcceleration");
  
  Param.GetScalar(IsothermalSoundSpeed, "Physics.Hydro.IsothermalSoundSpeed");
  
  Param.GetArray(ConstantAcceleration, "Physics.Gravity.ConstantAcceleration");
  Param.GetScalar(Mu, "Physics.Hydro.Mu");
  Param.GetScalar(DivBDampingLength, "Physics.MHD.DivBDampingLength");
  Param.GetScalar(CoolingCutOffDensity1, "Physics.AtomicPhysics.CoolingCutOffDensity1");
  Param.GetScalar(CoolingCutOffDensity2, "Physics.AtomicPhysics.CoolingCutOffDensity2");
  Param.GetScalar(CoolingCutOffTemperature, "Physics.AtomicPhysics.CoolingCutOffTemperature");
  Param.GetScalar(CoolingPowerCutOffDensity1, "Physics.AtomicPhysics.CoolingPowerCutOffDensity1");
  Param.GetScalar(CoolingPowerCutOffDensity2, "Physics.AtomicPhysics.CoolingPowerCutOffDensity2");
  Param.GetScalar(UseH2OnDust, "Physics.AtomicPhysics.UseH2OnDust");
  Param.GetScalar(UseCUDA, "SimulationControl.UseCUDA");
  
  Param.GetScalar(MoveParticlesBetweenSiblings, "SimulationControl.Optimization.MoveParticlesBetweenSiblings");
  Param.GetScalar(ParticleSplitterIterations, "SimulationControl.Optimization.ParticleSplitterIterations");
  Param.GetScalar(ParticleSplitterChildrenParticleSeparation, "SimulationControl.Optimization.ParticleSplitterChildrenParticleSeparation");
  Param.GetScalar(ResetMagneticField, "Physics.MHD.ResetMagneticField");
  Param.GetArray(ResetMagneticFieldAmplitude, "Physics.MHD.ResetMagneticFieldAmplitude");

  
  /* Now we know which hydro solver we're using, we can assign the
     default Riemann solver and flux reconstruction methods.  These
     parameters aren't used for PPM_LagrangeRemap and Zeus. */
  
  if (HydroMethod == PPM_DirectEuler) {
    if (RiemannSolver == INT_UNDEFINED) 
      RiemannSolver = TwoShock;
    if (ReconstructionMethod == INT_UNDEFINED)
      ReconstructionMethod = PPM;
    if (RiemannSolver == HLL && ReconstructionMethod == PLM) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("RiemannSolver = HLL, ReconstructionMethod = PLM.\n"
	       "These are the defaults for the MUSCL (hydro_rk) solvers,\n"
	       "but don't exist for the FORTRAN solvers.  Changing to the\n"
	       "old FORTRAN default, two-shock and PPM.\n"
	       "-- To override this, set RiemannSolver = -1\n");
      RiemannSolver = TwoShock;
      ReconstructionMethod = PPM;
    }
    if (RiemannSolver == -HLL) RiemannSolver = HLL;
  } else if (HydroMethod == HD_RK || HydroMethod == MHD_RK) {
    if (RiemannSolver == INT_UNDEFINED) 
      RiemannSolver = HLL;
    if (ReconstructionMethod == INT_UNDEFINED)
      ReconstructionMethod = PLM;
  }
  
  //  OutputTemperature = ((ProblemType == 7) || (ProblemType == 11));
  
  /* Even if this is not cosmology, due to a check for nested grid cosmology
     in ProtoSubgrid_AcceptableGrid.C, we'll set the default for this here. */
  CosmologySimulationNumberOfInitialGrids = 1;
  
  if (HydroMethod != MHD_RK)
    BAnyl = 0; // set this to zero no matter what unless we have a magnetic field to analyze.
  
  /* Count static nested grids since this isn't written in the
     parameter file */
  
  // for (i = 0; i < MAX_STATIC_REGIONS; i++)
  //   if (StaticRefineRegionLevel[i] != INT_UNDEFINED)
  //     CosmologySimulationNumberOfInitialGrids++;

  CosmologySimulationNumberOfInitialGrids += NumberOfStaticRefineRegions;
	
	
  /* If we have turned on Comoving coordinates, read cosmology parameters. */
	
  if (ComovingCoordinates) {
    
    // Always output temperature in cosmology runs
    OutputTemperature = TRUE;
    
    if (CosmologyReadParameters(&MetaData.StopTime, &MetaData.Time)
	== FAIL) {
      ENZO_FAIL("Error in ReadCosmologyParameters.\n");
    }
  }
  else {
    if (ReadUnits() == FAIL){
      ENZO_FAIL("Error in ReadUnits. ");
    }
  }
	
  // make sure that MHD is turned on if we're trying to use anisotropic conduction.
  // if not, alert user.
  if(AnisotropicConduction==TRUE && useMHD==0){
    ENZO_FAIL("AnisotropicConduction can only be used if MHD is turned on!\n");
  }  
  if(AnisotropicConduction==TRUE && MetaData.TopGridRank < 2){
    ENZO_FAIL("AnisotropicConduction can only be used if TopGridRank is >= 2!\n");
  }
  

  /*
    if (EOSType == 3) // an isothermal equation of state implies the adiabatic index = 1 
    Gamma = 1; 
  */

  /* convert MustRefineParticlesMinimumMass from a mass into a density, 
     ASSUMING CUBE simulation space */
  MustRefineParticlesMinimumMass /= POW(1/(float(MetaData.TopGridDims[0])*POW(float(RefineBy), float(MustRefineParticlesRefineToLevel))),3);
  
  
  /* Use Physical units stuff */
  
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0, PressureUnits = 1.0;
  double MassUnits = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, 
	     &MassUnits, MetaData.Time);
    PressureUnits = DensityUnits*pow(VelocityUnits,2);
    /*IMOPORTANT: If change anything here must change both equivilant parts in WriteParameterFile.C as well */
    
    /* Change input physical parameters into code units */
    MustRefineParticlesMinimumMass /= MassUnits; 
    StarMakerOverDensityThreshold /= DensityUnits;
    //  StarEnergyFeedbackRate = StarEnergyFeedbackRate/pow(LengthUnits,2)*pow(TimeUnits,3);
    
    if (SinkMergeDistance > 1.0)
      SinkMergeDistance /= LengthUnits;
    //printf(" \n SinkMergeDistance = %"FSYM"\n \n", SinkMergeDistance);
    SmallRho /= DensityUnits;
    SmallP /= PressureUnits;
    SmallT /= TemperatureUnits;
    MaximumAlvenSpeed /= VelocityUnits;
    EOSSoundSpeed /=  VelocityUnits;
    float h, cs, dpdrho, dpde;
    EOS(SmallP, SmallRho, SmallEint, h, cs, dpdrho, dpde, EOSType, 1);
    if (debug && (HydroMethod == HD_RK || HydroMethod == MHD_RK))
      printf("smallrho=%g, smallp=%g, smalleint=%g, DensityUnits = %g, PressureUnits=%g, MaximumAlvenSpeed=%g\n",
	     SmallRho, SmallP, SmallEint,DensityUnits, PressureUnits, MaximumAlvenSpeed);
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) 
      if (MinimumMassForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumMassForRefinement[i] /= MassUnits;
      }
    if (GravitationalConstant > 12.49 && GravitationalConstant < 12.61) {
      GravitationalConstant = 4.0 * 3.1415926 * 6.6726e-8 * DensityUnits * pow(TimeUnits,2);
      printf("Gravitational Constant recalculated from 4pi to 4piG in code units\n");
    }
    
  }
  
  /* For !restart, this only ruins the units because
     MinimumOverDensityForRefinement is already set in defaults.cfg
     and not FLOAT_UNDEFINED.  For restart,
     MinimumOverDensityForRefinement is not even needs to be read
     because only MinimumMassForRefinement is used for CellFlagging.
     So, why did we have to do this in the first place?  - Ji-hoon Kim
     in Apr.2010 (The counterpart in WriteParameterFile is also
     commented out) */   //#####
  
  /*
  if (!ComovingCoordinates && UsePhysicalUnit) 
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) 
      if (MinimumOverDensityForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumOverDensityForRefinement[i] /= DensityUnits;
      }
  */
  
  /* If RefineRegionTimeType is 0 or 1, read in the input file. */
  if ((RefineRegionTimeType == 0) || (RefineRegionTimeType == 1)) {
    if (ReadEvolveRefineFile() == FAIL) {
      ENZO_FAIL("Error in ReadEvolveRefineFile.");
    }
  }
  
  /* If GadgetEquilibriumCooling == TRUE, we don't want MultiSpecies
     or RadiationFieldType to be on - both are taken care of in
     the Gadget cooling routine.  Therefore, we turn them off!
     Also, initialize the Gadget equilibrium cooling data. */
  
  if(GadgetEquilibriumCooling == TRUE){
    
    if(MyProcessorNumber == ROOT_PROCESSOR ) {
      fprintf(stderr, "WARNING:  GadgetEquilibriumCooling = 1.  Forcing\n");
      fprintf(stderr, "WARNING:  RadiationFieldType = 0, MultiSpecies = 0, and\n");
      fprintf(stderr, "WARNING:  RadiativeCooling = 1.\n");
    }
    
    RadiationFieldType = 0;
    MultiSpecies       = 0;
    RadiativeCooling   = 1;

    // initialize Gadget equilibrium cooling
    if (InitializeGadgetEquilibriumCoolData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeGadgetEquilibriumCoolData.");
    } 
  }
  
  /* If set, initialize the RadiativeCooling and RateEquations data. */
  
  if (MultiSpecies > 0) {
    if (InitializeRateData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeRateData.");
    }
  }
  
  if (MultiSpecies             == 0 && 
      MetalCooling             == 0 &&
      GadgetEquilibriumCooling == 0 &&
      RadiativeCooling          > 0) {
    if (InitializeEquilibriumCoolData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeEquilibriumCoolData.");
    }
  }
  
  /* If set, initialze Cosmic Ray Efficiency Models */
  
  if (CRModel){
    if (InitializeCosmicRayData() == FAIL){
      ENZO_FAIL("Error in Initialize CosmicRayData.");
    }
  }
  
  /* If using the internal radiation field, initialize it. */
  
  if (RadiationFieldType == 11) 
    RadiationData.RadiationShield = TRUE; 
  else if (RadiationFieldType == 10)
    RadiationData.RadiationShield = FALSE; 
  
  if ((RadiationFieldType >= 10 && RadiationFieldType <= 11) ||
      RadiationData.RadiationShield == TRUE)
    if (InitializeRadiationFieldData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeRadiationFieldData.");
    }
  
  /* If using MBHFeedback = 2 to 5 (Star->FeedbackFlag = MBH_JETS), 
     you need MBHParticleIO for angular momentum */
  
  if (MBHFeedback >= 2 && MBHFeedback <= 5) 
    MBHParticleIO = TRUE;
  
  /* Turn off DualEnergyFormalism for zeus hydro (and a few other things). */
  
  if (HydroMethod == Zeus_Hydro) {
    ConservativeInterpolation = FALSE;
    DualEnergyFormalism       = FALSE;
    //    FluxCorrection            = FALSE;
  }
  
  /* Set some star feedback parameters. */
  
  if ((STARFEED_METHOD(NORMAL_STAR) || STARFEED_METHOD(UNIGRID_STAR)) && 
      (StarFeedbackDistRadius > 0)) {
    
    // Calculate number of cells in the shape over which to distribute feedback.
    StarFeedbackDistRadius = min(StarFeedbackDistRadius,
				 StarFeedbackDistCellStep);
    int i, j, k, cell_step;
    
    StarFeedbackDistTotalCells = 0;
    for (k = -StarFeedbackDistRadius;k <= StarFeedbackDistRadius;k++) {
      for (j = -StarFeedbackDistRadius;j <= StarFeedbackDistRadius;j++) {
	for (i = -StarFeedbackDistRadius;i <= StarFeedbackDistRadius;i++) {
	  cell_step = fabs(k) + fabs(j) + fabs(i);
	  if (cell_step <= StarFeedbackDistCellStep) {
	    StarFeedbackDistTotalCells++;
	  }
	}
      }
    }
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr,"Total cells for star feedback smoothing: %"ISYM".\n",
	      StarFeedbackDistTotalCells);
    }
  }
  
  /* For rk_hydro, we need to set some variables */
  
  if (DualEnergyFormalism) {
    NEQ_HYDRO = 6;
    NEQ_MHD   = 10;
    ieint = 5;
    iBx = 6;
    iBy = 7;
    iBz = 8;
    iPhi = 9;
    iEint = 5;
  }
  
  // Don't include free electron field
  switch (MultiSpecies) {
  case 0:  NSpecies = 0;  break;
  case 1:  NSpecies = 5;  break;
  case 2:  NSpecies = 8;  break;
  case 3:  NSpecies = 11; break;
  default: NSpecies = 0;  break;
  }
  
  // Determine color fields (NColor) later inside a grid object.
  // ...
  
  /* Set the number of particle attributes, if left unset. */
  
  if (NumberOfParticleAttributes == INT_UNDEFINED)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfParticleAttributes = 3;
    else
      NumberOfParticleAttributes = 0;
  
#ifdef UNUSED
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = (RadiativeCooling && SelfGravity
				     && HydroMethod == Zeus_Hydro) ?
      max(MaximumRefinementLevel-2, 5) : MaximumRefinementLevel;
#else
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = MaximumRefinementLevel;
#endif
  
  
  MaximumGravityRefinementLevel =
    min(MaximumGravityRefinementLevel, MaximumRefinementLevel);
  
  /* If MultiSpecies < 2, we can't simulate Pop III star formation */
  
  if (MultiSpecies < 2 && STARMAKE_METHOD(POP3_STAR)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Cannot form Pop III stars without H2 cooling!\n"
	      "Turning Pop III star formation OFF.\n");
    StarParticleCreation -= 1 << POP3_STAR;
  }
  
  /* Use the value in MaximumParticleRefinementLevel to set the smoothing
     radius for the particles, to be used to Grid_DepositPositions. */
  
  if (MaximumParticleRefinementLevel >= 0)
    DepositPositionsParticleSmoothRadius =
      (DomainRightEdge[0] - DomainLeftEdge[0])/
      (float(MetaData.TopGridDims[0])*
       POW(float(RefineBy), float(MaximumParticleRefinementLevel)));
  else
    DepositPositionsParticleSmoothRadius = 0;
  
  
  
  //  PPMDiffusion causes an out-of-bounds condition as currently written
  //  The following is an over-ride to force PPMDiffusion OFF. This has
  //  been fixed in this latest version (AK).
  
  //  NOTE: The fix keeps the code from crashing, but is not a proper 
  //  implementation of PPM diffusion.  The reason why is that Enzo typically
  //  uses 3 ghost zones, and the correct PPM diffusion implementation requires
  //  4 parameters.  SO, you should not use this parameter for, e.g., cosmology
  //  runs unless you know what you're doing.  (BWO)
  
  if (MetaData.PPMDiffusionParameter != 0 && ProblemType != 60 // Turbulence
      && ProblemType != 4  // Double Mach Reflection test
      && ProblemType != 6  // Implosion test
      && ProblemType != 7  // SedovBlast test
      && ProblemType != 8  // KH test
      && ProblemType != 9  // Noh test
      && ProblemType != 11 // Radiating shock test
      ) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("WARNING! Setting MetaData.PPMDiffusionParameter = 0\n");
    MetaData.PPMDiffusionParameter = 0;
  }
  
  if (PartitionNestedGrids == 1) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("WARNING! PartitionNestedGrids = 1 forces Parallel IO = 1\n");
    ParallelRootGridIO = 1;
    ParallelParticleIO = 1;
  }
  
  if(ProblemType==70 && UseHydro==1){
    printf("ReadParameterFile: ProblemType=70.  Disabling hydrodynamics!\n");
    UseHydro=FALSE;
  }
  
  
  if ((MetaData.GravityBoundary != TopGridPeriodic) &&
      (UnigridTranspose)) {
    /* it turns out that Robert Harkness' unigrid transpose stuff is incompatible with the top
       grid isolated gravity boundary conditions.  I'm not 100 percent sure why this is - in the 
       meantime, just double-check to make sure that if one tries to use the isolated boundary
       conditions when the unigrid transpose stuff is on, the code crashes loudly.
       -- BWO, 26 June 2008 */
    if (MyProcessorNumber == ROOT_PROCESSOR){
      fprintf(stderr, "\n\n");
      fprintf(stderr, "  ************************************************************************\n");
      fprintf(stderr, "  ****  D'oh!  At present, you cannot use isolated top grid boundary  ****\n");
      fprintf(stderr, "  ****  conditions with the top grid unigrid bookkeeping scheme.      ****\n");
      fprintf(stderr, "  ****  Consult Brian O'Shea for the details of this wackiness,       ****\n");
      fprintf(stderr, "  ****  and in the meantime enzo DISABLED unigrid tranposition!       ****\n");
      fprintf(stderr, "  ************************************************************************\n");      
      fprintf(stderr, "\n\n");
    }
    UnigridTranspose = FALSE;
  }
  
  /* If the restart dump parameters were set to the previous defaults
     (dtRestartDump = 5 hours), then set back to current default,
     which is no restart dumps. */
  
  float tol = 1e-6;
  if (ABS(MetaData.dtRestartDump - 3600.0*5) / (3600.0*5) < tol &&
      ABS(MetaData.TimeLastRestartDump) < tol) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	      "==================================================\n"
	      "-> Turning off restart dumps because the previous\n"
	      "default was set to 5 hours but not used.  To avoid this warning,\n"
	      "set dtRestartDump to a negative value.  If you wanted to use\n"
	      "restart dumps, please set it != 18000.\n"
	      "==================================================\n");
    MetaData.dtRestartDump = FLOAT_UNDEFINED;
  }
  
  /* If refining by must-refine particles, particle mass refinement
     must be turned on. */
  
  int method;
  bool TurnOnParticleMassRefinement = false;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) 
    if (CellFlaggingMethod[method] == 8) {
      TurnOnParticleMassRefinement = true;
      break;
    }
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) 
    if (CellFlaggingMethod[method] == 4) {
      TurnOnParticleMassRefinement = false;
      break;
    }
  if (TurnOnParticleMassRefinement) {
    method = 0;
    while (CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }
  
  
  /* If we're refining the region around P3 supernovae,
     MustRefineByParticles must be set.  Check this.  */
  
  if (PopIIISupernovaMustRefine == TRUE) {
    bool TurnOnParticleMustRefine = true;
    for (method = 0; method < MAX_FLAGGING_METHODS; method++)
      if (CellFlaggingMethod[method] == 8)
	TurnOnParticleMustRefine = false;
    if (TurnOnParticleMustRefine) {
      method = 0;
      while (CellFlaggingMethod[method] != INT_UNDEFINED)
	method++;
      CellFlaggingMethod[method] = 8;
    }
    
    /* Check if the must refine level is still at the default.  If so,
       break because it's zero!  Won't do anything, and the user will
       be disappointed to find that the simulation didn't refine
       around the SN. */
    
    if (MustRefineParticlesRefineToLevel == 0)
      ENZO_FAIL("MustRefineParticlesRefineToLevel is still ZERO, and you set"
		"PopIIISupernovaMustRefine.  Set the level or turn off"
		"PopIIISupernovaMustRefine.");
  } // ENDIF PopIIISupernovaMustRefine
  
  if (TracerParticleOn) {
    ParticleTypeInFile = TRUE;
  }
  
  if (OutputParticleTypeGrouping && (!ParticleTypeInFile)) {
    OutputParticleTypeGrouping = FALSE;
  }
  
  //  if (WritePotential && ComovingCoordinates && SelfGravity) {
  if (WritePotential && SelfGravity) {
    CopyGravPotential = TRUE;
  }
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    
    if ( MetaData.GlobalDir != NULL ) {
      fprintf(stderr, "Output to Global Dir %s\n", MetaData.GlobalDir);
    }
    
    if ( MetaData.LocalDir != NULL ) {
      fprintf(stderr, "Output to Local Dir %s\n", MetaData.LocalDir);
    }
    
  }
  
  // mqk temporary:
  MetaData.LocalDir = NULL;
  if ( (MetaData.GlobalDir != NULL) && (MetaData.LocalDir != NULL) ) {
    ENZO_FAIL("Cannot select GlobalDir AND LocalDir!\n");
  }
  
  char *cwd_buffer = new char[MAX_LINE_LENGTH];
  size_t cwd_buffer_len = MAX_LINE_LENGTH;
  
  if ( (MetaData.GlobalDir == NULL) && (MetaData.LocalDir == NULL) ) {
    /*if(getcwd(cwd_buffer, cwd_buffer_len) == NULL) {
      fprintf(stderr, "GETCWD call FAILED\n");
      }
      if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"CWD %s\n", cwd_buffer);
    */
    /* No one seems to want GlobalDir to default to abspath(CWD).  I'm leaving
       the code here in case you do. MJT */ 
    strcpy(cwd_buffer, ".");
    MetaData.GlobalDir = cwd_buffer;
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"Global Dir set to %s\n", cwd_buffer);
  }
  
  /* Generate unique identifier if one wasn't found. */
  if(MetaData.SimulationUUID == NULL){
    MetaData.SimulationUUID = new char[MAX_LINE_LENGTH];
    get_uuid(MetaData.SimulationUUID);
  }
  
  for (int i=0; i<MetaData.TopGridRank;i++)
    TopGridDx[i]=(DomainRightEdge[i]-DomainLeftEdge[i])/MetaData.TopGridDims[i];
  
  //  for (int i=0; i<MetaData.TopGridRank; i++)
  //      fprintf (stderr, "read  %"ISYM"  %"ISYM" \n", 
  // 	      MetaData.LeftFaceBoundaryCondition[i], 
  // 	      MetaData.RightFaceBoundaryCondition[i]);
  
  if (UseCUDA) {
    LoadBalancing = 0; // Should explore how LoadBalancing = 1 gives problems with CUDA
#ifndef ECUDA
    printf("This executable was compiled without CUDA support.\n");
    printf("use \n");
    printf("make cuda-yes\n");
    printf("Exiting.\n");
    my_exit(EXIT_SUCCESS);
#endif
  }
  if (debug) printf("Initialdt in ReadParameterFile = %e\n", *Initialdt);
  
  CheckShearingBoundaryConsistency(MetaData);
	
	
  delete dummy;
	
  return SUCCESS;
#endif /* ndef CONFIG_USE_LIBCONFIG */
}
