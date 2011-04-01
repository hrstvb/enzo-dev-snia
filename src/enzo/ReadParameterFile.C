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

#include "parameters.h"


/* This variable is declared here and only used in Grid_ReadGrid. */
 


/* function prototypes */

void my_exit(int status); 
int ReadListOfFloats(FILE *fptr, int N, float floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime);
int ReadUnits(FILE *fptr);
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
 
  /* read MetaData parameters */
 
  Param.GetScalar(MetaData.CycleNumber, "InternalParameters.InitialCycleNumber");
  Param.GetScalar(MetaData.Time, "InternalParameters.InitialTime");
  Param.GetScalar(MetaData.CPUTime, "InternalParameters.InitialCPUTime");
  Param.GetScalar((*Initialdt), "Initialdt");
  
  Param.GetScalar(CheckpointRestart, "CheckpointRestart"); // should be bool
  Param.GetScalar(MetaData.StopTime, "SimulationControl.StopTime");
  Param.GetScalar(MetaData.StopCycle, "SimulationControl.StopCycle");
  Param.GetScalar(MetaData.StopSteps, "SimulationControl.StopSteps");
  Param.GetScalar(MetaData.StopCPUTime, "StopCPUTime");
  Param.GetScalar(MetaData.ResubmitOn, "SimulationControl.ResubmitOn"); // should be bool
  
  
  Param.GetScalar(MetaData.ResubmitCommand, "SimulationControl.ResubmitCommand");
  
  Param.GetScalar(MetaData.MaximumTopGridTimeStep, "SimulationControl.MaximumTopGridTimeStep");
  
  Param.GetScalar(MetaData.TimeLastRestartDump, "InternalParameters.outputLabeling.TimeLastRestartDump");
  Param.GetScalar(MetaData.dtRestartDump, "OutputControlParameters.restartDump.dtRestartDump");
  Param.GetScalar(MetaData.TimeLastDataDump, "InternalParameters.outputLabeling.TimeLastDataDump");
  Param.GetScalar(MetaData.dtDataDump, "OutputControlParameters.dataDump.dtDataDump");
  Param.GetScalar(MetaData.TimeLastHistoryDump, "InternalParameters.outputLabeling.TimeLastHistoryDump");
  Param.GetScalar(MetaData.dtHistoryDump, "dtHistoryDump");
  
  Param.GetScalar(TracerParticleOn, "SimulationControl.TracerParticleOn"); // should be bool
  Param.GetScalar(ParticleTypeInFile, "OutputControlParameters.ParticleTypeInFile"); // should be bool
  Param.GetScalar(OutputParticleTypeGrouping, "OutputControlParameters.OutputParticleTypeGrouping");
  Param.GetScalar(MetaData.TimeLastTracerParticleDump, "TimeLastTracerParticleDump"); //doc?
  Param.GetScalar(MetaData.dtTracerParticleDump, "dtTracerParticleDump");//doc?
  Param.GetScalar(MetaData.TimeLastInterpolatedDataDump, "TimeLastInterpolatedDataDump");//doc?
  Param.GetScalar(MetaData.dtInterpolatedDataDump, "dtInterpolatedDataDump");//doc?
  
  
  Param.GetArray(MetaData.NewMovieLeftEdge, "OutputControlParameters.movieDump.NewMovieLeftEdge");
  Param.GetArray(MetaData.NewMovieRightEdge, "OutputControlParameters.movieDump.NewMovieLeftEdge");
  
  Param.GetScalar(MetaData.CycleLastRestartDump, "CycleLastRestartDump"); //not used
  Param.GetScalar(MetaData.CycleSkipRestartDump, "CycleSkipRestartDump"); //not used
  Param.GetScalar(MetaData.CycleLastDataDump, "OutputControlParameters.outputLabeling.CycleLastDataDump");
  Param.GetScalar(MetaData.CycleSkipDataDump, "OutputControlParameters.cycleDump.CycleSkipDataDump");
  Param.GetScalar(MetaData.CycleLastHistoryDump, "OutputControlParameters.outputLabeling.CycleLastHistoryDump");
  Param.GetScalar(MetaData.CycleSkipHistoryDump, "OutputControlParameters.cycleDump.CycleSkipHistoryDump");
  Param.GetScalar(MetaData.CycleSkipGlobalDataDump, "CycleSkipGlobalDataDump");
  Param.GetScalar(MetaData.OutputFirstTimeAtLevel, "OutputControlParameters.outputTriggers.OutputFirstTimeAtLevel");
  Param.GetScalar(MetaData.StopFirstTimeAtLevel, "SimulationControl.StopFirstTimeAtLevel");
  
  
  /* Maximum density directed output */
  Param.GetScalar(OutputOnDensity, "OutputControlParameters.outputTriggers.OutputOnDensity");  // should be bool
  Param.GetScalar(StartDensityOutputs, "OutputControlParameters.outputTriggers.StartDensityOutputs");
  Param.GetScalar(CurrentDensityOutput, "OutputControlParameters.outputTriggers.CurrentDensityOutput");
  IncrementDensityOutput = Param.GetScalar("OutputControlParameters.outputTriggers.IncrementDensityOutput");
  
  
  /* Subcycle directed output */
  Param.GetScalar(MetaData.SubcycleSkipDataDump, "OutputControlParameters.cycleDump.SubcycleSkipDataDump");
  Param.GetScalar(MetaData.SubcycleLastDataDump, "InternalParameters.outputLabeling.SubcycleLastDataDump");
  Param.GetScalar(MetaData.SubcycleNumber, "InternalParameters.SubcycleNumber");
  
  
  Param.GetScalar(FileDirectedOutput, "OutputControlParameters.FileDirectedOutput");
  Param.GetScalar(WriteBinaryHierarchy, "OutputControlParameters.WriteBinaryHierarchy");
  
  
  Param.GetScalar(MetaData.RestartDumpNumber, "InternalParameters.outputLabeling.RestartDumpNumber");
  Param.GetScalar(MetaData.DataDumpNumber, "InternalParameters.outputLabeling.DataDumpNumber");
  Param.GetScalar(MetaData.HistoryDumpNumber, "InternalParameters.outputLabeling.HistoryDumpNumber");
  Param.GetScalar(MetaData.TracerParticleDumpNumber, "InternalParameters.outputLabeling.TracerParticleDumpNumber");
  
  Param.GetScalar(MetaData.RestartDumpName, "OutputControlParameters.restartDump.RestartDumpName");
  Param.GetScalar(MetaData.RestartDumpDir, "OutputControlParameters.restartDump.RestartDumpDir");
  
  Param.GetScalar(MetaData.DataDumpName, "OutputControlParameters.dataDump.DataDumpName");
  Param.GetScalar(MetaData.DataDumpDir, "OutputControlParameters.dataDump.DataDumpDir");
  
  Param.GetScalar(MetaData.RedshiftDumpName, "OutputControlParameters.restartDump.RedshiftDumpName");
  Param.GetScalar(MetaData.RedshiftDumpDir, "OutputControlParameters.redshiftDump.RedshiftDumpDir");
  
  Param.GetScalar(MetaData.TracerParticleDumpName, "OutputControlParameters.tracerParticleDump.TracerParticleDumpName");
  Param.GetScalar(MetaData.TracerParticleDumpDir, "OutputControlParameters.tracerParticleDump.TracerParticleDumpDir");
  
  Param.GetScalar(MetaData.HistoryDumpName, "OutputControlParameters.historyDump.HistoryDumpName");
  Param.GetScalar(MetaData.HistoryDumpDir, "OutputControlParameters.historyDump.HistoryDumpDir");
 
  Param.GetScalar(MetaData.LocalDir, "LocalDir");
  Param.GetScalar(MetaData.GlobalDir, "GlobalDir");
  
  Param.GetScalar(LoadBalancing, "SimulationControl.optimization.LoadBalancing");
  Param.GetScalar(ResetLoadBalancing, "SimulationControl.optimization.ResetLoadBalancing");
  Param.GetScalar(LoadBalancingCycleSkip, "SimulationControl.optimization.LoadBalancingCycleSkip");
  Param.GetScalar(LoadBalancingMinLevel, "SimulationControl.optimization.LoadBalancingMinLevel");
  Param.GetScalar(LoadBalancingMaxLevel, "SimulationControl.optimization.LoadBalancingMaxLevel");
  
  
  int NumberOfTimeActions = Param.Size("SimulationControl.timeaction.actions");
  if (NumberOfTimeActions > MAX_TIME_ACTIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of TimeActions (%d)!\n",MAX_TIME_ACTIONS)
      }
  
  char TimeActionNames[MAX_LINE_LENGTH][NumberOfTimeActions];
  Param.GetArray(TimeActionNames, "SimulationControl.timeaction.actions");
  
  for (i = 0; i < NumberOfTimeActions; i++) {
    Param.GetScalar(TimeActionType[i], "SimulationControl.timeaction.%s.Type", TimeActionNames[i]);
    Param.GetScalar(TimeActionRedshift[i], "SimulationControl.timeaction.%s.Redshift", TimeActionNames[i]);
    Param.GetScalar(TimeActionTime[i], "SimulationControl.timeaction.%s.Time", TimeActionNames[i]);
    Param.GetScalar(TimeActionParameter[i], "SimulationControl.timeaction.%s.Parameter", TimeActionParameter[i]);
  }
  
  Param.GetScalar(MetaData.StaticHierarchy, "SimulationControl.amr.StaticHierarchy");  // should be bool
  
  Param.GetScalar(MetaData.TopGridRank, "SimulationControl.domain.TopGridRank");
  Param.GetArray(MetaData.TopGridDims, "TopGridDimensions");
  
  Param.GetScalar(MetaData.GravityBoundary, "PhysicsParameters.TopGridGravityBoundary");
  
#ifdef TRANSFER
  Param.GetScalar(MetaData.RadHydroParameterFname, "RadHydroParamfile");
#endif
  Param.GetScalar(ImplicitProblem, "PhysicsParameters.radiationTransfer.ImplicitProblem"); // should be bool
  Param.GetScalar(RadiativeTransferFLD, "PhysicsParameters.radiationTransfer.RadiativeTransferFLD");
#ifdef EMISSIVITY
  Param.GetScalar(StarMakerEmissivityField, "StarMakerEmissivityField");
  Param.GetScalar(uv_param, "uv_param");
#endif
  
  Param.GetScalar(MetaData.ParticleBoundaryType, "PhysicsParameters.ParticleBoundaryType");
  Param.GetScalar(MetaData.NumberOfParticles, "InternalParameters.NumberOfParticles");
  
  Param.GetScalar(MetaData.CourantSafetyNumber, "PhysicsParameters.hydro.CourantSafetyNumber");
  Param.GetScalar(MetaData.PPMFlatteningParameter, "PhysicsParameters.hydro.PPMFlatteningParameter");  // should be bool
  Param.GetScalar(MetaData.PPMDiffusionParameter, "PhysicsParameters.hydro.PPMDiffusionParameter");  // should be bool
  Param.GetScalar(MetaData.PPMSteepeningParameter, "PhysicsParameters.hydro.PPMSteepeningParameter");  // should be bool
  
  
  /* read global Parameters */
  Param.GetScalar(ProblemType, "Initialization.ProblemType");
  
#ifdef NEW_PROBLEM_TYPES
  Param.GetScalar(ProblemTypeName, "Initialization.ProblemTypeName");
  ProblemType = -978;
#endif
  
  Param.GetScalar(HydroMethod, "PhysicsParameters.hydro.HydroMethod");
  
  if (HydroMethod==MHD_RK) useMHD = 1;
  
  
  Param.GetScalar(huge_number, "SimulationControl.huge_number");
  Param.GetScalar(tiny_number, "SimulationControl.tiny_number");
  param.GetScalar(Gamma, "PhysicsParameters.hydro.Gamma");
  Param.GetScalar(PressureFree, "PhysicsParameters.hydro.PressureFree");
  Param.GetScalar(RefineBy, "SimulationControl.amr.RefineBy");
  Param.GetScalar(MaximumRefinementLevel, "SimulationControl.amr.MaximumRefinementLevel");
  Param.GetScalar(MaximumGravityRefinementLevel, "SimulationControl.amr.MaximumGravityRefinementLevel");
  Param.GetScalar(MaximumParticleRefinementLevel, "SimulationControl.amr.MaximumParticleRefinementLevel");
  Param.GetArray("SimulationControl.amr.CellFlaggingMethod", CellFlaggingMethod);
  
  Param.GetScalar(FluxCorrection, "PhysicsParameters.hydro.Fluxcorrection");
  Param.GetScalar(InterpolationMethod, "PhysicsParameters.hydro.InterpolationMethod");
  Param.GetScalar(ConservativeInterpolation, "PhysicsParameters.hydro.ConservativeInterpolation");
  Param.GetScalar(MinimumEfficiency, "PhysicsParameters.optimization.MinimumEfficiency");
  Param.GetScalar(SubgridSizeAutoAdjust, "PhysicsParameters.optimization.SubgridSizeAutoAdjust");
  Param.GetScalar(OptimalSubgridsPerProcessor, "PhysicsParameters.optimization.OptimalSubgridsPerProcessor");
  Param.GetScalar(MinimumSubgridEdge, "PhysicsParameters.optimization.MinimumSubgridEdge");
  Param.GetScalar(MaximumSubgridSize, "PhysicsParameters.optimization.MaximumSubgridSize");
  Param.GetScalar(NumberOfBufferZones, "PhysicsParameters.amr.NumberOfBufferZones");
  Param.GetScalar(FastSiblingLocatorEntireDomain, "PhysicsParameters.optimization.FastSiblingLocatorEntireDomain");
  Param.GetScalar(MustRefineRegionMinRefinementLevel, "PhysicsParameters.amr.MustRefineRegionMinRefinementLevel");
  Param.GetScalar(MetallicityRefinementMinLevel, "PhysicsParameters.amr.MetallicityRefinementMinLevel");  
  Param.GetScalar(MetallicityRefinementMinMetallicity, "PhysicsParameters.amr.MetallicityRefinementMinMetallicity"); 
  Param.GetScalar(MetallicityRefinementMinDensity, "PhysicsParameters.amr.MetallicityRefinementMinDensity"); 
  Param.GetArray(DomainLeftEdge, "SimulationControl.domain.DomainLeftEdge");
  Param.GetArray(DomainRightEdge, "SimulationControl.domain.DomainRightEdge");
  Param.GetArray(GridVelocity, "PhysicsParameters.GridVelocity");
  
  Param.GetScalar(RefineRegionAutoAdjust, "SimulationControl.domain.RefineRegionAutoAdjust"); 
  Param.GetArray(RefineRegionLeftEdge, "SimulationControl.domain.RefineRegionLeftEdge");
  Param.GetArray(RefineRegionRightEdge, "SimulationControl.domain.RefineRegionRightEdge");
  Param.GetArray(MustRefineRegionLeftEdge, "SimulationControl.domain.MustRefineRegionLeftEdge");
  Param.GetArray(MustRefineRegionRightEdge, "SimulationControl.domain.MustRefineRegionRightEdge");
  
  
  /* Read evolving RefineRegion */
  
  Param.GetScalar(RefineRegionTimeType, "SimulationControl.amr.RefineRegionTimeType");
  Param.GetScalar(RefineRegionFile, "RefineRegionFile");
  
  int NumberOfFields = Param.Size("InternalParameters.fields");
  if (NumberOfFields > MAX_NUMBER_OF_BARYON_FIELDS) {
    ENZO_VFAIL("You've exceeded the maximum number of BaryonFields (%d)!\n",MAX_NUMBER_OF_BARYON_FIELDS)
      }
  
  char FieldNames[MAX_LINE_LENGTH][NumberOfFields];
  Param.GetArray(StaticRefineRegionNames, "InternalParameters.fields");
  
  for (i = 0; i < NumberOfFields; i++) {
    Param.GetScalar(DataLabel[i], "InternalParameters.%s.name", FieldNames[i]);
    Param.GetScalar(DataUnits[i], "InternalParameters.%s.cgsConversionFactor", FieldNames[i]);
  }
  
  Param.GetScalar(UniformGravity, "PhysicsParameters.gravity.UniformGravity");
  Param.GetScalar(UniformGravityDirection, "PhysicsParameters.gravity.UniformGravityDirection");
  Param.GetScalar(UniformGravityConstant, "PhysicsParameters.gravity.UniformGravityConstant");
  
  Param.GetScalar(PointSourceGravity, "PhysicsParameters.gravity.PointSourceGravity"); // should be bool
  Param.GetArray(PointSourceGravityPosition, "PhysicsParameters.gravity.PointSourceGravityPosition");
  Param.GetScalar(PointSourceGravityConstant, "PhysicsParameters.gravity.PointSourceGravityConstant");
  Param.GetScalar(PointSourceGravityCoreRadius, "PhysicsParameters.gravity.PointSourceGravityCoreRadius");
  
  Param.GetScalar(ExternalGravity, "PhysicsParameters.gravity.ExternalGravity");
  
  Param.GetScalar(SelfGravity, "PhysicsParameters.gravity.SelfGravity");
  Param.GetScalar(SelfGravityGasOff, "PhysicsParameters.gravity.SelfGravityGasOff");
  Param.GetScalar(AccretionKernel, "PhysicsParameters.AccretionKernel");
  Param.GetScalar(GravitationalConstant, "PhysicsParameters.gravity.GravitationalConstant");
  
  Param.GetScalar(S2ParticleSize, "S2ParticleSize");
  Param.GetScalar(GravityResolution, "PhysicsParameters.gravity.GravityResolution");
  Param.GetScalar(ComputePotential, "PhysicsParameters.gravity.ComputePotential"); // should be bool
  Param.GetScalar(PotentialIterations, "PhysicsParameters.gravity.PotentialIterations");
  Param.GetScalar(WritePotential, "OutputControlParameters.supplementalFields.WritePotential");  // should be bool
  Param.GetScalar(BaryonSelfGravityApproximation, "PhysicsParameters.gravity.BaryonSelfGravityApproximation");  // should be bool
  
  Param.GetScalar(GreensFunctionMaxNumber, "PhysicsParameters.gravity.GreensFunctionMaxNumber");
  Param.GetScalar(GreensFunctionMaxSize, "PhysicsParameters.gravity.GreensFunctionMaxSize");
  
  Param.GetScalar(DualEnergyFormalism, "PhysicsParameters.hydro.DualEnergyFormalism"); // should be bool
  Param.GetScalar(DualEnergyFormalismEta1, "PhysicsParameters.hydro.DualEnergyFormalismEta1");
  Param.GetScalar(DualEnergyFormalismEta2, "PhysicsParameters.hydro.DualEnergyFormalismEta2");
  Param.GetScalar(ParticleCourantSafetyNumber, "PhysicsParameters.hydro.ParticleCourantSafetyNumber");
  Param.GetScalar(RootGridCourantSafetyNumber, "PhysicsParameters.hydro.RootGridCourantSafetyNumber");
  Param.GetScalar(RandomForcing, "PhysicsParameters.miscellaneous.RandomForcing");
  Param.GetScalar(RandomForcingEdot, "PhysicsParameters.miscellaneous.RandomForcingEdot");
  Param.GetScalar(RandomForcingMachNumber, "RandomForcingMachNumber");
  Param.GetScalar(RadiativeCooling, "RadiativeCooling"); // should be bool
  Param.GetScalar(RadiativeCoolingModel, "PhysicsParameters.atomicPhysics.RadiativeCoolingModel");
  Param.GetScalar(GadgetEquilibriumCooling, "PhysicsParameters.atomicPhysics.GadgetEquilibriumCooling");
  Param.GetScalar(MultiSpecies, "PhysicsParameters.atomicPhysics.MultiSpecies");
  Param.GetScalar(PrimordialChemistrySolver, "PhysicsParameters.atomicPhysics.PrimordialChemistrySolver");
  Param.GetScalar(CIECooling, "PhysicsParameters.atomicPhysics.CIECooling"); // should be bool
  Param.GetScalar(H2OpticalDepthApproximation, "PhysicsParameters.atomicPhysics.H2OpticalDepthApproximation"); // should be bool
  Param.GetScalar(ThreeBodyRate, "PhysicsParameters.atomicPhysics.ThreeBodyRate"); // should be bool
  
  
  Param.GetScalar(CloudyCoolingData.CloudyCoolingGridFile, "PhysicsParameters.atomicPhysics.cloudyCooling");
  Param.GetScalar(CloudyCoolingData.IncludeCloudyHeating, "PhysicsParameters.atomicPhysics.cloudyCooling.IncludeCloudyHeating"); // should be bool
  Param.GetScalar(CloudyCoolingData.IncludeCloudyMMW, "PhysicsParameters.atomicPhysics.cloudyCooling.IncludeCloudyMMW"); // should be bool
  Param.GetScalar(CloudyCoolingData.CMBTemperatureFloor, "PhysicsParameters.atomicPhysics.cloudyCooling.CMBTemperatureFloor"); // should be bool
  Param.GetScalar(CloudyCoolingData.CloudyMetallicityNormalization, "PhysicsParameters.atomicPhysics.cloudyCooling.CloudyMetallicityNormalization");
  Param.GetScalar(CloudyCoolingData.CloudyElectronFractionFactor, "PhysicsParameters.atomicPhysics.cloudyCooling.CloudyElectronFractionFactor");
  Param.GetScalar(MetalCooling, "PhysicsParameters.atomicPhysics.MetalCooling"); // should be bool
  
  Param.GetScalar(MetalCoolingTable, "PhysicsParameters.atomicPhysics.MetalCoolingTable");
  
  Param.GetScalar(CRModel, "PhysicsParameters.miscellaneous.CRModel");
  Param.GetScalar(ShockMethod, "PhysicsParameters.miscellaneous.ShockMethod");
  Param.GetScalar(ShockTemperatureFloor, "PhysicsParameters.miscellaneous.ShockTemperatureFloor");
  Param.GetScalar(StorePreShockFields, "PhysicsParameters.miscellaneous.StorePreShockFields");
  
  Param.GetScalar(RadiationFieldType, "PhysicsParameters.radiationField.RadiationFieldType");
  Param.GetScalar(TabulatedLWBackground, "TabulatedLWBackground"); // should be bool
  Param.GetScalar(AdjustUVBackground, "PhysicsParameters.radiationField.AdjustUVBackground");
  Param.GetScalar(SetUVBAmplitude, "PhysicsParameters.radiationField.SetUVBAmplitude");
  Param.GetScalar(SetHeIIHeatingScale, "PhysicsParameters.radiationField.SetHeIIHeatingScale");
  Param.GetScalar(RadiationFieldLevelRecompute, "PhysicsParameters.radiationField.RadiationFieldLevelRecompute"); // should be bool
  Param.GetScalar(RadiationData.RadiationShield, "PhysicsParameters.radiationField.RadiationShield"); // should be bool
  Param.GetScalar(CoolData.f3, "PhysicsParameters.radiationField.RadiationSpectrumNormalization");
  Param.GetScalar(CoolData.alpha0, "PhysicsParameters.radiationField.RadiationSpectrumSlope");
  Param.GetScalar(CoolData.f0to3, "PhysicsParameters.atomicPhysics.CoolDataf0to3");
  Param.GetScalar(CoolData.RadiationRedshiftOn, "PhysicsParameters.radiationField.RadiationRedshiftOn");
  Param.GetScalar(CoolData.RadiationRedshiftOff, "PhysicsParameters.radiationField.RadiationRedshiftOff");
  Param.GetScalar(CoolData.RadiationRedshiftFullOn, "PhysicsParameters.radiationField.RadiationRedshiftFullOn");
  Param.GetScalar(CoolData.RadiationRedshiftDropOff, "PhysicsParameters.radiationField.RadiationRedshiftDropOff");
  Param.GetScalar(CoolData.HydrogenFractionByMass, "PhysicsParameters.atomicPhysics.HydrogenFractionByMass");
  Param.GetScalar(CoolData.DeuteriumToHydrogenRatio, "PhysicsParameters.atomicPhysics.DeuteriumToHydrogenRatio");
  Param.GetScalar(CoolData.NumberOfTemperatureBins, "PhysicsParameters.atomicPhysics.NumberOfTemperatureBins");
  Param.GetScalar(CoolData.ih2co, "PhysicsParameters.atomicPhysics.CoolDataIh2co"); // should be bool
  Param.GetScalar(CoolData.ipihdt, "PhysicsParameters.atomicPhysics.CoolDataIpiht"); // should be bool
  Param.GetScalar(CoolData.TemperatureStart, "PhysicsParameters.atomicPhysics.TemperatureStart");
  Param.GetScalar(CoolData.TemperatureEnd, "PhysicsParameters.atomicPhysics.TemperatureEnd");
  Param.GetScalar(CoolData.comp_xray, "PhysicsParameters.atomicPhysics.CoolDataCompXray"); // should be bool
  Param.GetScalar(CoolData.temp_xray, "PhysicsParameters.atomicPhysics.CoolDataTempXray"); // should be bool
  Param.GetScalar(RateData.CaseBRecombination, "RateDataCaseBRecombination"); // should be bool
  Param.GetScalar(PhotoelectricHeating, "PhysicsParameters.atomicPhysics.PhotoelectricHeating"); // should be bool
  
  Param.GetScalar(OutputCoolingTime, "OutputControlParameters.supplementalFields.OutputCoolingTime"); // should be bool 
  Param.GetScalar(OutputTemperature, "OutputControlParameters.supplementalFields.OutputTemperature"); // should be bool 
  Param.GetScalar(OutputSmoothedDarkMatter, "OutputControlParameters.supplementalFields.OutputSmoothedDarkMatter"); // should be bool 
  Param.GetScalar(SmoothedDarkMatterNeighbors, "OutputControlParameters.supplementalFields.SmoothedDarkMatterNeighbors"); // should be bool 
  Param.GetScalar(OutputGriddedStarParticle, "OutputControlParameters.supplementalFields.OutputGriddedStarParticle"); // should be bool 
  
  Param.GetScalar(ZEUSQuadraticArtificialViscosity, "PhysicsParameters.hydro.ZEUSQuadraticArtificialViscosity");
  Param.GetScalar(ZEUSLinearArtificialViscosity, "PhysicsParameters.hydro.ZEUSLinearArtificialViscosity");
  
  Param.GetScalar(UseMinimumPressureSupport, "PhysicsParameters.hydro.UseMinimumPressureSupport"); // should be bool 
  Param.GetScalar(MinimumPressureSupportParameter, "PhysicsParameters.hydro.MinimumPressureSupportParameter");
  Param.GetScalar(RefineByJeansLengthSafetyFactor, "SimulationControl.amr.RefineByJeansLengthSafetyFactor");
  Param.GetScalar(RefineByJeansLengthUnits, "SimulationControl.amr.RefineByJeansLengthUnits"); // should be bool 
  Param.GetScalar(JeansRefinementColdTemperature, "SimulationControl.amr.JeansRefinementColdTemperature");
  Param.GetScalar(RefineByResistiveLengthSafetyFactor, "SimulationControl.amr.RefineByResistiveLengthSafetyFactor");
  Param.GetScalar(MustRefineParticlesRefineToLevel, "SimulationControl.amr.MustRefineParticlesRefineToLevel");
  Param.GetScalar(MustRefineParticlesRefineToLevelAutoAdjust, "SimulationControl.amr.MustRefineParticlesRefineToLevelAutoAdjust"); // should be bool 
  Param.GetScalar(MustRefineParticlesMinimumMass, "SimulationControl.amr.MustRefineParticlesMinimumMass");
  Param.GetScalar(ParticleTypeInFile, "OutputControlParameters.ParticleTypeInFile"); // should be bool 
  
  int NumberOfStaticRefineRegions = Param.Size("StaticRefineRegion.regions");
  if (NumberOfStaticRefineRegions > MAX_STATIC_REGIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of StaticRefineRegions (%d)!\n",MAX_STATIC_REGIONS)
      }
  
  char StaticRefineRegionNames[MAX_LINE_LENGTH][NumberOfStaticRefineRegions];
  Param.GetArray(StaticRefineRegionNames, "StaticRefineRegion.regions");
  
  for (i = 0; i < NumberOfStaticRefineRegions; i++) {
    Param.GetScalar(StaticRefineRegionLevel[i], "StaticRefineRegion.%s.Level",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionLeftEdge[i],"StaticRefineRegion.%s.LeftEdge",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionRightEdge[i],"StaticRefineRegion.%s.RightEdge",StaticRefineRegionNames[i]);
  }
  
  Param.GetScalar(ParallelRootGridIO, "InternalParameters.outputLabeling.ParallelRootGridIO"); // should be bool
  Param.GetScalar(ParallelParticleIO, "InternalParameters.outputLabeling.ParallelParticleIO"); // should be bool
  
  Param.GetScalar(Unigrid, "InternalParameters.outputLabeling.Unigrid"); // should be bool
  Param.GetScalar(UnigridTranspose, "InternalParameters.outputLabeling.UnigridTranspose"); // should be bool
  Param.GetScalar(NumberOfRootGridTilesPerDimensionPerProcessor, "InternalParameters.outputLabeling.NumberOfRootGridTilesPerDimensionPerProcessor");
  
  Param.GetScalar(PartitionNestedGrids, "Initialization.PartitionNestedGrids");
  Param.GetScalar(StaticPartitionNestedGrids, "Initialization.StaticPartitionNestedGrids"); // should be bool
  
  Param.GetScalar(ExtractFieldsOnly, "OutputControlParameters.ExtractFieldsOnly"); // should be bool
  
  Param.GetScalar(debug1, "InternalParameters.Debug1"); // should be bool
  
  Param.GetScalar(debug2, "InternalParameters.Debug2"); // should be bool
  
  Param.GetScalar(MemoryLimit, "SimulationControl.optimization.MemoryLimit");
  
#ifdef STAGE_INPUT
  Param.GetScalar(StageInput, "StageInput");
#endif
  
#ifdef OOC_BOUNDARY
  
  Param.GetScalar(ExternalBoundaryIO, "ExternalBoundaryIO"); // should be bool
  
  Param.GetScalar(ExternalBoundaryTypeIO, "ExternalBoundaryTypeIO"); // should be bool
  
  Param.GetScalar(ExternalBoundaryValueIO, "ExternalBoundaryValueIO"); // should be bool
  
  Param.GetScalar(SimpleConstantBoundary, "SimpleConstantBoundary"); // should be bool
  
#endif
  
  Param.GetArray(SlopeFlaggingFields, "SimulationControl.amr.SlopeFlaggingFields");
  
  Param.GetArray(MinimumSlopeForRefinement, "SimulationControl.amr.MinimumSlopeForRefinement");
  
  Param.GetArray(MinimumOverDensityForRefinement, "SimulationControl.amr.MinimumOverDensityForRefinement");
  
  Param.GetArray(MinimumMassForRefinement, "SimulationControl.amr.MinimumMassForRefinement");
  
  Param.GetArray(MinimumMassForRefinementLevelExponent, "SimulationControl.amr.MinimumMassForRefinementLevelExponent");
  
  
  Param.GetScalar(MinimumPressureJumpForRefinement, "SimulationControl.amr.MinimumPressureJumpForRefinement");
  Param.GetScalar(MinimumShearForRefinement, "SimulationControl.amr.MinimumShearForRefinement");
  Param.GetScalar(MinimumEnergyRatioForRefinement, "SimulationControl.amr.MinimumEnergyRatioForRefinement");
  Param.GetScalar(ShockwaveRefinementMinMach, "SimulationControl.amr.ShockwaveRefinementMinMach");
  Param.GetScalar(ShockwaveRefinementMinVelocity, "SimulationControl.amr.ShockwaveRefinementMinVelocity");
  Param.GetScalar(ShockwaveRefinementMaxLevel, "SimulationControl.amr.ShockwaveRefinementMaxLevel");
  Param.GetScalar(ComovingCoordinates, "PhysicsParameters.cosmology.ComovingCoordinates"); // should be bool
  Param.GetScalar(StarParticleCreation, "PhysicsParameters.otherParticles.StarParticleCreation");
  Param.GetScalar(BigStarFormation, "PhysicsParameters.otherParticles.bigStar.BigStarFormation"); // should be bool
  Param.GetScalar(BigStarFormationDone, "PhysicsParameters.otherParticles.bigStar.BigStarFormationDone"); // should be bool
  Param.GetScalar(BigStarSeparation, "PhysicsParameters.otherParticles.bigStar.BigStarSeparation");
  Param.GetScalar(SimpleQ, "SimulationControl.radiationTransfer.SimpleQ");
  Param.GetScalar(SimpleRampTime, "SimulationControl.radiationTransfer.SimpleRampTime");
  Param.GetScalar(StarParticleFeedback, "PhysicsParameters.otherParticles.StarParticleFeedback");
  Param.GetScalar(NumberOfParticleAttributes, "PhysicsParameters.otherParticles.NumberOfParticleAttributes");
  Param.GetScalar(AddParticleAttributes, "PhysicsParameters.otherParticles.AddParticleAttributes"); // should be bool
  
  
  /* read data which defines the boundary conditions */
  Param.GetArray(MetaData.LeftFaceBoundaryCondition, "PhysicsParameters.LeftFaceBoundaryCondition");
  Param.GetArray(MetaData.RightFaceBoundaryCondition, "PhysicsParameters.RightFaceBoundaryCondition");
  
  Param.GetScalar(MetaData.BoundaryConditionName, "BoundaryConditionName");
  Param.GetScalar(MetaData.MetaDataIdentifier, "MetaDataIdentifier");
  Param.GetScalar(MetaData.SimulationUUID, "MetaDataSimulationUUID");
  Param.GetScalar(MetaData.RestartDatasetUUID, "MetaDataDatasetUUID");
  Param.GetScalar(MetaData.InitialConditionsUUID, "MetaDataInitialConditionsUUID");
  
  /* Check version number. */
  
  Param.GetScalar(TempFloat, "InternalParameters.provenance.VersionNumber");
  if (fabs(TempFloat - VERSION) >= 1.0e-3 &&
      MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Warning: Incorrect version number.\n");
  
  
  /* Read star particle parameters. */
  
  Param.GetScalar(StarMakerOverDensityThreshold, "PhysicsParameters.otherParticles.starMaker.StarMakerOverDensityThreshold");
  Param.GetScalar(StarMakerSHDensityThreshold, "PhysicsParameters.otherParticles.starMaker.StarMakerSHDensityThreshold");
  Param.GetScalar(StarMakerMassEfficiency, "PhysicsParameters.otherParticles.starMaker.StarMakerMassEfficiency");
  Param.GetScalar(StarMakerMinimumMass, "PhysicsParameters.otherParticles.starMaker.StarMakerMinimumMass");
  Param.GetScalar(StarMakerMinimumDynamicalTime, "PhysicsParameters.otherParticles.starMaker.StarMakerMinimumDynamicalTime");
  Param.GetScalar(StarMassEjectionFraction, "PhysicsParameters.otherParticles.starMaker.StarMassEjectionFraction");
  Param.GetScalar(StarMetalYield, "PhysicsParameters.otherParticles.starMaker.StarMetalYield");
  Param.GetScalar(StarEnergyToThermalFeedback, "PhysicsParameters.otherParticles.starEnergy.StarEnergyToThermalFeedback");
  Param.GetScalar(StarEnergyToStellarUV, "PhysicsParameters.otherParticles.starEnergy.StarEnergy.ToStellarUV");
  Param.GetScalar(StarEnergyToQuasarUV, "PhysicsParameters.otherParticles.starEnergy.StarEnergyToQuasarUV");
  
  Param.GetScalar(StarFeedbackDistRadius, "StarFeedbackDistRadius");
  Param.GetScalar(StarFeedbackDistCellStep, "StarFeedbackDistCellStep");
  Param.GetScalar(StarClusterUseMetalField, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterUseMetalField");
  Param.GetScalar(StarClusterMinDynamicalTime, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterMinDynamicalTime");
  Param.GetScalar(StarClusterHeliumIonization, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterHeliumIonization");
  Param.GetScalar(StarClusterUnresolvedModel, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterUnresolvedModel");
  Param.GetScalar(StarClusterIonizingLuminosity, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterIonizingLuminosity");
  Param.GetScalar(StarClusterSNEnergy, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterSNEnergy");
  Param.GetScalar(StarClusterSNRadius, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterSNRadius");
  Param.GetScalar(StarClusterFormEfficiency, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterFormEfficiency");
  Param.GetScalar(StarClusterMinimumMass, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterMinimumMass");
  Param.GetScalar(StarClusterCombineRadius, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterCombineRadius");
  
  Param.GetArray(StarClusterRegionLeftEdge, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterRegionLeftEdge");
  Param.GetArray(StarClusterRegionRightEdge, "PhysicsParameters.otherParticles.starClusterParticle.StarClusterRegionRightEdge");
  
  Param.GetScalar(PopIIIStarMass, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIStarMass");                           
  Param.GetScalar(PopIIIInitialMassFunctionSeed, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIInitialMassFunctionSeed");                 
  Param.GetScalar(PopIIIInitialMassFunctionCalls, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIInitialMassFunctionCalls"); 
  Param.GetScalar(PopIIILowerMassCutoff, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIILowerMassCutoff");
  Param.GetScalar(PopIIIUpperMassCutoff, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIUpperMassCutoff");
  
  //   ret += sscanf(line,"PopIIIMassRange %"FSYM" %"FSYM, &PopIIILowerMassCutoff, &PopIIIUpperMassCutoff);
  
  Param.GetScalar	(PopIIIInitialMassFunctionSlope, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIInitialMassFunctionSlope");
  Param.GetScalar            (PopIIIBlackHoles, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIBlackHoles");
  Param.GetScalar		 (PopIIIBHLuminosityEfficiency, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIBHLuminosityEfficiency");
  Param.GetScalar <float>		 (PopIIIOverDensityThreshold, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIOverDensityThreshold");
  Param.GetScalar <float>		(PopIIIH2CriticalFraction, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIH2CriticalFraction");
  Param.GetScalar <float>		(PopIIIMetalCriticalFraction, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIMetalCriticalFraction");
  Param.GetScalar          (PopIIISupernovaRadius, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIISupernovaRadius");
  Param.GetScalar		(PopIIISupernovaUseColour, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIISupernovaUseColour");
  Param.GetScalar <int>		(PopIIISupernovaMustRefine, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIISupernovaMustRefine");
  Param.GetScalar <int>		 (PopIIISupernovaMustRefineResolution, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIISupernovaMustRefineResolution");
  Param.GetScalar		 (PopIIIHeliumIonization, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIHeliumIonization");
  Param.GetScalar <float>		 (PopIIIColorDensityThreshold, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIColorDensityThreshold");
  Param.GetScalar <float>		 (PopIIIColorMass, "PhysicsParameters.otherParticles.popIIIStarParticle.PopIIIColorMass");
  
  Param.GetScalar(MBHAccretion, "PhysicsParameters.otherParticles.mbhParticle.MBHAccretion");
  Param.GetScalar(MBHAccretionRadius, "PhysicsParameters.otherParticles.mbhParticle.MBHAccretionRadius");
  Param.GetScalar(MBHAccretingMassRatio, "PhysicsParameters.otherParticles.mbhParticle.MBHAccretingMassRatio");
  Param.GetScalar(MBHAccretionFixedTemperature, "PhysicsParameters.otherParticles.mbhParticle.MBHAccretionFixedTemperature");
  Param.GetScalar(MBHAccretionFixedRate, "PhysicsParameters.otherParticles.mbhParticle.MBHAccretionFixedRate");
  Param.GetScalar(MBHTurnOffStarFormation, "PhysicsParameters.otherParticles.mbhParticle.MBHTurnOffStarFormation");
  Param.GetScalar(MBHCombineRadius, "PhysicsParameters.otherParticles.mbhParticle.MBHCombineRadius");
  Param.GetScalar(MBHMinDynamicalTime, "PhysicsParameters.otherParticles.mbhParticle.MBHMinDynamicalTime");
  Param.GetScalar(MBHMinimumMass, "PhysicsParameters.otherParticles.mbhParticle.MBHMinimumMass");
  
  Param.GetScalar(MBHFeedback, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedback");
  Param.GetScalar(MBHFeedbackRadiativeEfficiency, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackRadiativeEfficiency");
  Param.GetScalar(MBHFeedbackEnergyCoupling, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackEnergyCoupling");
  Param.GetScalar(MBHFeedbackMassEjectionFraction, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackMassEjectionFraction");
  Param.GetScalar(MBHFeedbackMetalYield, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackMetalYield");
  Param.GetScalar(MBHFeedbackThermalRadius, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackThermalRadius");
  Param.GetScalar(MBHFeedbackJetsThresholdMass, "PhysicsParameters.otherParticles.mbhParticle.MBHFeedbackJetsThresholdMass");
  Param.GetScalar <int>		  (MBHParticleIO, "PhysicsParameters.otherParticles.mbhParticle.MBHParticleIO");
  
  Param.GetScalar(MBHParticleIOFilename, "MBHParticleIOFilename");
  Param.GetScalar(MBHInsertLocationFilename, "MBHInsertLocationFilename");
  
  
  /* Read Movie Dump parameters */
  
  Param.GetScalar(MovieSkipTimestep, "OutputControlParameters.movieDump.MovieSkipTimestep");
  Param.GetScalar(Movie3DVolumes, "OutputControlParameters.movieDump.Movie3DVolumes");
  Param.GetScalar(MovieVertexCentered, "OutputControlParameters.movieDump.MovieVertexCentered");
  Param.GetScalar(NewMovieParticleOn, "OutputControlParameters.movieDump.NewMovieParticleOn");
  Param.GetArray(MovieDataField, "MovieDataField");
  Param.GetScalar(NewMovieDumpNumber, "OutputControlParameters.movieDump.NewMovieDumpNumber");
  Param.GetScalar(NewMovieName, "NewMovieName");
  Param.GetScalar(MetaData.MovieTimestepCounter, "OutputControlParameters.movieDump.MovieTimestepCounter");
  
  Param.GetScalar(MultiMetals, "PhysicsParameters.atomicPhysics.MultiMetals");
  Param.GetScalar(IsotropicConduction, "IsotropicConduction"); // should be bool
  Param.GetScalar(AnisotropicConduction, "AnisotropicConduction"); // should be bool
  Param.GetScalar(IsotropicConductionSpitzerFraction, "IsotropicConductionSpitzerFraction");
  Param.GetScalar(AnisotropicConductionSpitzerFraction, "AnisotropicConductionSpitzerFraction");
  Param.GetScalar(ConductionCourantSafetyNumber, "PhysicsParameters.conduction.ConductionCourantSafetyNumber");
  
  Param.GetScalar(RadiativeTransfer, "PhysicsParameters.radiationTransfer.RadiativeTransfer"); // should be bool
  Param.GetScalar(RadiationXRaySecondaryIon, "PhysicsParameters.radiationTransfer.RadiationXRaySecondaryIon"); // should be bool
  Param.GetScalar(RadiationXRayComptonHeating, "PhysicsParameters.radiationTransfer.RadiationXRayComptonHeating"); // should be bool
  
  
  /* Shearing Box Boundary parameters */
  
  Param.GetScalar(AngularVelocity, "PhysicsParameters.AngularVelocity");
  Param.GetScalar(VelocityGradient, "PhysicsParameters.VelocityGradient");
  Param.GetScalar(ShearingVelocityDirection, "PhysicsParameters.ShearingVelocityDirection");
  Param.GetScalar(ShearingBoxProblemType, "PhysicsParameters.ShearingBoxProblemType");
  
  
// #ifdef STAGE_INPUT
//   sscanf(line, "LocalPath = %s\n", LocalPath);
//   sscanf(line, "GlobalPath = %s\n", GlobalPath);
// #endif
  
  /* Embedded Python */
  Param.GetScalar(PythonTopGridSkip, "AnalysisParameters.python.PythonTopGridSkip"); // should be bool
  Param.GetScalar(PythonSubcycleSkip, "AnalysisParameters.python.PythonSubcycleSkip"); // should be bool
  
#ifdef USE_PYTHON
  Param.GetScalar(NumberOfPythonCalls, "NumberOfPythonCalls");
  Param.GetScalar(NumberOfPythonTopGridCalls, "NumberOfPythonTopGridCalls");
  Param.GetScalar(NumberOfPythonSubcycleCalls, "NumberOfPythonSubcycleCalls");
#endif
  
  /* Inline halo finder */
  
  Param.GetScalar(InlineHaloFinder, "AnalysisParameters.haloFinder.InlineHaloFinder"); // should be bool
  Param.GetScalar(HaloFinderSubfind, "AnalysisParameters.haloFinder.HaloFinderSubfind"); // should be bool
  Param.GetScalar(HaloFinderOutputParticleList, "AnalysisParameters.haloFinder.HaloFinderOutputParticleList"); // should be bool
  Param.GetScalar(HaloFinderRunAfterOutput, "AnalysisParameters.haloFinder.HaloFinderRunAfterOutput"); // should be bool
  Param.GetScalar(HaloFinderLinkingLength, "AnalysisParameters.haloFinder.HaloFinderLinkingLength");
  Param.GetScalar(HaloFinderMinimumSize, "AnalysisParameters.haloFinder.HaloFinderMinimumSize");
  Param.GetScalar(HaloFinderCycleSkip, "AnalysisParameters.haloFinder.HaloFinderCycleSkip");
  Param.GetScalar(HaloFinderTimestep, "AnalysisParameters.haloFinder.HaloFinderTimestep");
  Param.GetScalar(HaloFinderLastTime, "AnalysisParameters.haloFinder.HaloFinderLastTime");
  
  /* This Block for Stanford Hydro */
  
  Param.GetScalar(UseHydro, "PhysicsParameters.hydro.UseHydro"); // should be bool
  
  
  /* Sink particles (for present day star formation) & winds */
  Param.GetScalar(SinkMergeDistance, "PhysicsParameters.otherParticles.sinkParticle.SinkMergeDistance");
  Param.GetScalar(SinkMergeMass, "PhysicsParameters.otherParticles.sinkParticle.SinkMergeMass");
  Param.GetScalar(StellarWindFeedback, "PhysicsParameters.otherParticles.StellarWindFeedback"); // should be bool
  Param.GetScalar(StellarWindTurnOnMass, "PhysicsParameters.otherParticles.StellarWindTurnOnMass");
  Param.GetScalar(MSStellarWindTurnOnMass, "PhysicsParameters.otherParticles.MSStellarWindTurnOnMass");
  
  Param.GetScalar(VelAnyl, "OutputControlParameters.supplementalFields.VelAnyl"); // should be bool
  Param.GetScalar(BAnyl, "OutputControlParameters.supplementalFields.BAnyl"); // should be bool
  
  
  
  /* Read MHD Paramters */
  Param.GetScalar(UseDivergenceCleaning, "PhysicsParameters.mhd.UseDivergenceCleaning"); // should be bool
  Param.GetScalar(DivergenceCleaningBoundaryBuffer, "PhysicsParameters.mhd.DivergenceCleaningBoundaryBuffer"); // should be bool
  Param.GetScalar(DivergenceCleaningThreshold, "PhysicsParameters.mhd.DivergenceCleaningThreshold");
  Param.GetScalar(PoissonApproximationThreshold, "PhysicsParameters.PoissonApproximationThreshold");
  Param.GetScalar(PoissonBoundaryType, "PhysicsParameters.PoissonBoundaryType");
  
  
  Param.GetScalar(AngularVelocity, "PhysicsParameters.AngularVelocity");
  Param.GetScalar(VelocityGradient, "PhysicsParameters.VelocityGradient");
  Param.GetScalar(UseDrivingField, "PhysicsParameters.miscellaneous.UseDrivingField"); // should be bool
  Param.GetScalar(DrivingEfficiency, "PhysicsParameters.miscellaneous.DrivingEfficiency");
  
  Param.GetScalar(StringKick, "StringKick");
  Param.GetScalar(StringKickDimension, "StringKickDimension");
  Param.GetScalar(UsePhysicalUnit, "PhysicsParameters.UsePhysicalUnit"); // should be bool
  Param.GetScalar(Theta_Limiter, "PhysicsParameters.hydro.Theta_Limiter");
  Param.GetScalar(RKOrder, "PhysicsParameters.hydro.RKOrder");
  Param.GetScalar(UseFloor, "PhysicsParameters.hydro.UseFloor"); // should be bool
  Param.GetScalar(UseViscosity, "PhysicsParameters.hydro.UseViscosity"); // should be bool
  Param.GetScalar(ViscosityCoefficient, "PhysicsParameters.hydro.ViscosityCoefficient");
  
  Param.GetScalar(UseAmbipolarDiffusion, "PhysicsParameters.hydro.UseAmbipolarDiffusion"); // should be bool
  Param.GetScalar(UseResistivity, "PhysicsParameters.hydro.UseResistivity"); // should be bool
  Param.GetScalar(SmallRho, "PhysicsParameters.SmallRho");
  Param.GetScalar(SmallP, "PhysicsParameters.SmallP");
  Param.GetScalar(SmallT, "PhysicsParameters.SmallT");
  Param.GetScalar(MaximumAlvenSpeed, "PhysicsParameters.mhd.MaximumAlvenSpeed");
  Param.GetScalar(Coordinate, "PhysicsParameters.Coordinate");
  Param.GetScalar(RiemannSolver, "PhysicsParameters.hydro.RiemannSolver");
  Param.GetScalar(RiemannSolverFallback, "PhysicsParameters.hydro.RiemannSolverFallback"); // should be bool
  Param.GetScalar(ConservativeReconstruction, "PhysicsParameters.hydro.ConservativeReconstruction"); // should be bool
  Param.GetScalar(PositiveReconstruction, "PhysicsParameters.hydro.PositiveReconstruction"); // should be bool
  Param.GetScalar(ReconstructionMethod, "PhysicsParameters.hydro.ReconstructionMethod");
  
  Param.GetScalar(EOSType, "PhysicsParameters.hydro.EOSType");
  Param.GetScalar(EOSSoundSpeed, "PhysicsParameters.hydro.EOSSoundSpeed");
  Param.GetScalar(EOSCriticalDensity, "PhysicsParameters.hydro.EOSCriticalDensity");
  Param.GetScalar(EOSGamma, "PhysicsParameters.hydro.EOSGamma");
  Param.GetScalar(UseConstantAcceleration, "PhysicsParameters.gravity.UseConstantAcceleration"); // should be bool
  
  Param.GetScalar(IsothermalSoundSpeed, "PhysicsParameters.hydro.IsothermalSoundSpeed");
  
  Param.GetArray(ConstantAcceleration, "PhysicsParameters.gravity.ConstantAcceleration");
  Param.GetScalar(Mu, "PhysicsParameters.hydro.Mu");
  Param.GetScalar(DivBDampingLength, "PhysicsParameters.mhd.DivBDampingLength");
  Param.GetScalar(CoolingCutOffDensity1, "PhysicsParameters.atomicPhysics.CoolingCutOffDensity1");
  Param.GetScalar(CoolingCutOffDensity2, "PhysicsParameters.atomicPhysics.CoolingCutOffDensity2");
  Param.GetScalar(CoolingCutOffTemperature, "PhysicsParameters.atomicPhysics.CoolingCutOffTemperature");
  Param.GetScalar(CoolingPowerCutOffDensity1, "PhysicsParameters.atomicPhysics.CoolingPowerCutOffDensity1");
  Param.GetScalar(CoolingPowerCutOffDensity2, "PhysicsParameters.atomicPhysics.CoolingPowerCutOffDensity2");
  Param.GetScalar(UseH2OnDust, "UseH2OnDust"); // should be bool
  Param.GetScalar(seCUDA, "UseCUDA"); // should be bool
  
  Param.GetScalar(MoveParticlesBetweenSiblings, "SimulationControl.optimization.MoveParticlesBetweenSiblings"); // should be bool
  Param.GetScalar(ParticleSplitterIterations, "SimulationControl.optimization.ParticleSplitterIterations");
  Param.GetScalar(ParticleSplitterChildrenParticleSeparation, "SimulationControl.optimization.ParticleSplitterChildrenParticleSeparation");
  Param.GetScalar(ResetMagneticField, "PhysicsParameters.mhd.ResetMagneticField"); // should be bool
  Param.GetArray(ResetMagneticFieldAmplitude, "PhysicsParameters.mhd.ResetMagneticFieldAmplitude");
  

#ifdef UNUSED
  /* check to see if the line belongs to one of the test problems */  

  if (strstr(line, "ShockTube")           ) ret++;
  if (strstr(line, "WavePool" )           ) ret++;
  if (strstr(line, "ShockPool")           ) ret++;
  if (strstr(line, "DoubleMach")          ) ret++;
  if (strstr(line, "Implosion")           ) ret++;
  if (strstr(line, "SedovBlast")          ) ret++;
  if (strstr(line, "Units")               ) ret++;
  if (strstr(line, "RadiatingShock")      ) ret++;
  if (strstr(line, "RotatingCylinder")    ) ret++;
  if (strstr(line, "RotatingSphere")    ) ret++;
  if (strstr(line, "StratifiedMediumExplosion")) ret++;
  if (strstr(line, "TestOrbit")    ) ret++;
  if (strstr(line, "KelvinHelmholtz")     ) ret++;
  if (strstr(line, "KH")                  ) ret++;
  if (strstr(line, "Noh")                 ) ret++;
  if (strstr(line, "TestProblem")         ) ret++;
  if (strstr(line, "ZeldovichPancake")    ) ret++;
  if (strstr(line, "PressurelessCollapse")) ret++;
  if (strstr(line, "AdiabaticExpansion" ) ) ret++;
  if (strstr(line, "CosmologySimulation") ) ret++;
  if (strstr(line, "TestGravity"        ) ) ret++;
  if (strstr(line, "SphericalInfall"    ) ) ret++;
  if (strstr(line, "TestGravitySphere"  ) ) ret++;
  if (strstr(line, "CollapseTest"       ) ) ret++;
  if (strstr(line, "Cosmology"          ) ) ret++;
  if (strstr(line, "SupernovaRestart"   ) ) ret++;
  if (strstr(line, "TracerParticleCreation")) ret++;
  if (strstr(line, "TurbulenceSimulation")) ret++;
  if (strstr(line, "ProtostellarCollapse")) ret++;
  if (strstr(line, "GalaxySimulation")) ret++;
  if (strstr(line, "ConductionTest")) ret++;
  if (strstr(line, "ConductionBubble")) ret++;
  if (strstr(line, "ConductionCloud")) ret++;
  if (strstr(line, "CoolingTest")) ret++;
  if (strstr(line, "ShearingBox")) ret++;
  if (strstr(line, "PoissonSolverTest")) ret++;
  /* 7.22.10 - CBH: Added 5 following lines to avoid runtime warnings from 
     extra params previously added to code (but not read_params) by others.*/
  if (strstr(line, "Cloudy")              ) ret++;
  if (strstr(line, "IsothermalSoundSpeed")) ret++;
  if (strstr(line, "dtPhoton")            ) ret++;
  if (strstr(line, "CurrentTimeIdentifier")) ret++;
  if (strstr(line, "MetaDataRestart")     ) ret++;
  if (strstr(line, "MustRefine") ) ret++;
  if (strstr(line, "AccretionKernel")     ) ret++;
  if (strstr(line, "PopIII")              ) ret++;
#ifdef TRANSFER
  if (strstr(line, "Radiative")           ) ret++;
  if (strstr(line, "PhotonTest")          ) ret++;
  
#endif
  
#endif  // #UNUSED


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

  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    if (StaticRefineRegionLevel[i] != INT_UNDEFINED)
      CosmologySimulationNumberOfInitialGrids++;
  
  /* If we have turned on Comoving coordinates, read cosmology parameters. */
 
  if (ComovingCoordinates) {
    
    // Always output temperature in cosmology runs
    OutputTemperature = TRUE;
    
    if (CosmologyReadParameters(fptr, &MetaData.StopTime, &MetaData.Time)
	== FAIL) {
      ENZO_FAIL("Error in ReadCosmologyParameters.\n");
    }
    rewind(fptr);
  }
  else {
    if (ReadUnits(fptr) == FAIL){
      ENZO_FAIL("Error in ReadUnits. ");
    }
    rewind(fptr);
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
  
  /* For !restart, this only ruins the units because MinimumOverDensityForRefinement is already 
     set in SetDefaultGlobalValues and not FLOAT_UNDEFINED.
     For restart, MinimumOverDensityForRefinement is not even needs to be read because only 
     MinimumMassForRefinement is used for CellFlagging.  
     So, why did we have to do this in the first place?  - Ji-hoon Kim in Apr.2010
     (The counterpart in WriteParameterFile is also commented out) */   //#####
  
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
  return SUCCESS;
#endif /* ndef CONFIG_USE_LIBCONFIG */
}
