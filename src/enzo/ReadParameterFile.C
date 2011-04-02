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
  
  Param.GetScalar(CheckpointRestart, "InternalParameters.outputLabeling.CheckpointRestart"); // should be bool
  Param.GetScalar(MetaData.StopTime, "SimulationControl.StopTime");
  Param.GetScalar(MetaData.StopCycle, "SimulationControl.StopCycle");
  Param.GetScalar(MetaData.StopSteps, "SimulationControl.StopSteps");
  Param.GetScalar(MetaData.StopCPUTime, "InternalParameters.StopCPUTime");
  Param.GetScalar(MetaData.ResubmitOn, "SimulationControl.ResubmitOn");
  
  Param.GetScalar(MetaData.ResubmitCommand, "SimulationControl.ResubmitCommand");
  
  Param.GetScalar(MetaData.MaximumTopGridTimeStep, "SimulationControl.MaximumTopGridTimeStep");
  
  Param.GetScalar(MetaData.TimeLastRestartDump, "InternalParameters.OutputLabeling.TimeLastRestartDump");
  Param.GetScalar(MetaData.dtRestartDump, "OutputControlParameters.RestartDump.dtRestartDump");
  Param.GetScalar(MetaData.TimeLastDataDump, "InternalParameters.OutputLabeling.TimeLastDataDump");
  Param.GetScalar(MetaData.dtDataDump, "OutputControlParameters.DataDump.dtDataDump");
  Param.GetScalar(MetaData.TimeLastHistoryDump, "InternalParameters.OutputLabeling.TimeLastHistoryDump");
  Param.GetScalar(MetaData.dtHistoryDump, "OutputControlParameters.HistoryDump.dtHistoryDump");
  
  Param.GetScalar(TracerParticleOn, "SimulationControl.TracerParticleOn");
  Param.GetScalar(ParticleTypeInFile, "OutputControlParameters.ParticleTypeInFile");
  Param.GetScalar(OutputParticleTypeGrouping, "OutputControlParameters.OutputParticleTypeGrouping");
  Param.GetScalar(MetaData.TimeLastTracerParticleDump, "InternalParameters.OutputLabeling.TimeLastTracerParticleDump"); //doc?
  Param.GetScalar(MetaData.dtTracerParticleDump, "OutputControlParameters.TracerParticleDump.dtTracerParticleDump");//doc?
  Param.GetScalar(MetaData.TimeLastInterpolatedDataDump, "InternalParameters.OutputLabeling.TimeLastInterpolatedDataDump");//doc?
  Param.GetScalar(MetaData.dtInterpolatedDataDump, "OutputControlParameters.dtInterpolatedDataDump");//doc?
  
  
  Param.GetArray(MetaData.NewMovieLeftEdge, "OutputControlParameters.MovieDump.NewMovieLeftEdge");
  Param.GetArray(MetaData.NewMovieRightEdge, "OutputControlParameters.MovieDump.NewMovieLeftEdge");
  
  Param.GetScalar(MetaData.CycleLastRestartDump, "OutputControlParameters.CycleDump.CycleLastRestartDump"); //not used
  Param.GetScalar(MetaData.CycleSkipRestartDump, "OutputControlParameters.CycleDump.CycleSkipRestartDump"); //not used
  Param.GetScalar(MetaData.CycleLastDataDump, "OutputControlParameters.OutputLabeling.CycleLastDataDump");
  Param.GetScalar(MetaData.CycleSkipDataDump, "OutputControlParameters.CycleDump.CycleSkipDataDump");
  Param.GetScalar(MetaData.CycleLastHistoryDump, "OutputControlParameters.OutputLabeling.CycleLastHistoryDump");
  Param.GetScalar(MetaData.CycleSkipHistoryDump, "OutputControlParameters.CycleDump.CycleSkipHistoryDump");
  Param.GetScalar(MetaData.CycleSkipGlobalDataDump, "AnalysisParameters.CycleSkipGlobalDataDump");
  Param.GetScalar(MetaData.OutputFirstTimeAtLevel, "OutputControlParameters.OutputTriggers.OutputFirstTimeAtLevel");
  Param.GetScalar(MetaData.StopFirstTimeAtLevel, "SimulationControl.StopFirstTimeAtLevel");
  
  
  /* Maximum density directed output */
  Param.GetScalar(OutputOnDensity, "OutputControlParameters.OutputTriggers.OutputOnDensity"); 
  Param.GetScalar(StartDensityOutputs, "OutputControlParameters.OutputTriggers.StartDensityOutputs");
  Param.GetScalar(CurrentDensityOutput, "OutputControlParameters.OutputTriggers.CurrentDensityOutput");
  IncrementDensityOutput = Param.GetScalar("OutputControlParameters.OutputTriggers.IncrementDensityOutput");
  
  
  /* Subcycle directed output */
  Param.GetScalar(MetaData.SubcycleSkipDataDump, "OutputControlParameters.CycleDump.SubcycleSkipDataDump");
  Param.GetScalar(MetaData.SubcycleLastDataDump, "InternalParameters.OutputLabeling.SubcycleLastDataDump");
  Param.GetScalar(MetaData.SubcycleNumber, "InternalParameters.SubcycleNumber");
  
  
  Param.GetScalar(FileDirectedOutput, "OutputControlParameters.FileDirectedOutput");
  Param.GetScalar(WriteBinaryHierarchy, "OutputControlParameters.WriteBinaryHierarchy");
  
  
  Param.GetScalar(MetaData.RestartDumpNumber, "InternalParameters.OutputLabeling.RestartDumpNumber");
  Param.GetScalar(MetaData.DataDumpNumber, "InternalParameters.OutputLabeling.DataDumpNumber");
  Param.GetScalar(MetaData.HistoryDumpNumber, "InternalParameters.OutputLabeling.HistoryDumpNumber");
  Param.GetScalar(MetaData.TracerParticleDumpNumber, "InternalParameters.OutputLabeling.TracerParticleDumpNumber");
  
  Param.GetScalar(MetaData.RestartDumpName, "OutputControlParameters.RestartDump.RestartDumpName");
  Param.GetScalar(MetaData.RestartDumpDir, "OutputControlParameters.RestartDump.RestartDumpDir");
  
  Param.GetScalar(MetaData.DataDumpName, "OutputControlParameters.DataDump.DataDumpName");
  Param.GetScalar(MetaData.DataDumpDir, "OutputControlParameters.DataDump.DataDumpDir");
  
  Param.GetScalar(MetaData.RedshiftDumpName, "OutputControlParameters.RestartDump.RedshiftDumpName");
  Param.GetScalar(MetaData.RedshiftDumpDir, "OutputControlParameters.RedshiftDump.RedshiftDumpDir");
  
  Param.GetScalar(MetaData.TracerParticleDumpName, "OutputControlParameters.TracerParticleDump.TracerParticleDumpName");
  Param.GetScalar(MetaData.TracerParticleDumpDir, "OutputControlParameters.TracerParticleDump.TracerParticleDumpDir");
  
  Param.GetScalar(MetaData.HistoryDumpName, "OutputControlParameters.HistoryDump.HistoryDumpName");
  Param.GetScalar(MetaData.HistoryDumpDir, "OutputControlParameters.HistoryDump.HistoryDumpDir");
  
  Param.GetScalar(MetaData.LocalDir, "OutputControlParameters.LocalDir");
  Param.GetScalar(MetaData.GlobalDir, "OutputControlParameters.GlobalDir");
  
  Param.GetScalar(LoadBalancing, "SimulationControl.Optimization.LoadBalancing");
  Param.GetScalar(ResetLoadBalancing, "SimulationControl.Optimization.ResetLoadBalancing");
  Param.GetScalar(LoadBalancingCycleSkip, "SimulationControl.Optimization.LoadBalancingCycleSkip");
  Param.GetScalar(LoadBalancingMinLevel, "SimulationControl.Optimization.LoadBalancingMinLevel");
  Param.GetScalar(LoadBalancingMaxLevel, "SimulationControl.Optimization.LoadBalancingMaxLevel");
  
  
  int NumberOfTimeActions = Param.Size("SimulationControl.Timeaction.Actions");
  if (NumberOfTimeActions > MAX_TIME_ACTIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of TimeActions (%d)!\n",MAX_TIME_ACTIONS)
      }
  
  char TimeActionNames[MAX_LINE_LENGTH][NumberOfTimeActions];
  Param.GetArray(TimeActionNames, "SimulationControl.Timeaction.Actions");
  
  for (i = 0; i < NumberOfTimeActions; i++) {
    Param.GetScalar(TimeActionType[i], "SimulationControl.Timeaction.%s.Type", TimeActionNames[i]);
    Param.GetScalar(TimeActionRedshift[i], "SimulationControl.Timeaction.%s.Redshift", TimeActionNames[i]);
    Param.GetScalar(TimeActionTime[i], "SimulationControl.Timeaction.%s.Time", TimeActionNames[i]);
    Param.GetScalar(TimeActionParameter[i], "SimulationControl.Timeaction.%s.Parameter", TimeActionParameter[i]);
  }
  
  Param.GetScalar(MetaData.StaticHierarchy, "SimulationControl.AMR.StaticHierarchy"); 
  
  Param.GetScalar(MetaData.TopGridRank, "SimulationControl.Domain.TopGridRank");
  Param.GetArray(MetaData.TopGridDims, "TopGridDimensions");
  
  Param.GetScalar(MetaData.GravityBoundary, "PhysicsParameters.TopGridGravityBoundary");
  
#ifdef TRANSFER
  Param.GetScalar(MetaData.RadHydroParameterFname, "PhysicsParameters.RadiationField.RadHydroParamfile");
#endif
  Param.GetScalar(ImplicitProblem, "PhysicsParameters.RadiationTransfer.ImplicitProblem");
  Param.GetScalar(RadiativeTransferFLD, "PhysicsParameters.RadiationTransfer.RadiativeTransferFLD");
#ifdef EMISSIVITY
  Param.GetScalar(StarMakerEmissivityField, "PhysicsParameters.OtherParticles.StarMaker.StarMakerEmissivityField");
  Param.GetScalar(uv_param, "PhysicsParameters.RadiationField.uv_param");
#endif
  
  Param.GetScalar(MetaData.ParticleBoundaryType, "PhysicsParameters.ParticleBoundaryType");
  Param.GetScalar(MetaData.NumberOfParticles, "InternalParameters.NumberOfParticles");
  
  Param.GetScalar(MetaData.CourantSafetyNumber, "PhysicsParameters.Hydro.CourantSafetyNumber");
  Param.GetScalar(MetaData.PPMFlatteningParameter, "PhysicsParameters.Hydro.PPMFlatteningParameter"); 
  Param.GetScalar(MetaData.PPMDiffusionParameter, "PhysicsParameters.Hydro.PPMDiffusionParameter"); 
  Param.GetScalar(MetaData.PPMSteepeningParameter, "PhysicsParameters.Hydro.PPMSteepeningParameter"); 
  
  
  /* read global Parameters */
  Param.GetScalar(ProblemType, "Initialization.ProblemType");
  
#ifdef NEW_PROBLEM_TYPES
  if (sscanf(line, "ProblemTypeName = %s", dummy) == 1) {
    ProblemTypeName = dummy;
    ProblemType = -978;
    ret = 1;
  }
#endif
  
  Param.GetScalar(HydroMethod, "PhysicsParameters.Hydro.HydroMethod");
  
  if (HydroMethod==MHD_RK) useMHD = 1;
  
  
  Param.GetScalar(huge_number, "SimulationControl.huge_number");
  Param.GetScalar(tiny_number, "SimulationControl.tiny_number");
  Param.GetScalar(Gamma, "PhysicsParameters.Hydro.Gamma");
  Param.GetScalar(PressureFree, "PhysicsParameters.Hydro.PressureFree");
  Param.GetScalar(RefineBy, "SimulationControl.AMR.RefineBy");
  Param.GetScalar(MaximumRefinementLevel, "SimulationControl.AMR.MaximumRefinementLevel");
  Param.GetScalar(MaximumGravityRefinementLevel, "SimulationControl.AMR.MaximumGravityRefinementLevel");
  Param.GetScalar(MaximumParticleRefinementLevel, "SimulationControl.AMR.MaximumParticleRefinementLevel");
  Param.GetArray("SimulationControl.AMR.CellFlaggingMethod", CellFlaggingMethod);
  
  Param.GetScalar(FluxCorrection, "PhysicsParameters.Hydro.Fluxcorrection");
  Param.GetScalar(InterpolationMethod, "PhysicsParameters.Hydro.InterpolationMethod");
  Param.GetScalar(ConservativeInterpolation, "PhysicsParameters.Hydro.ConservativeInterpolation");
  Param.GetScalar(MinimumEfficiency, "PhysicsParameters.Optimization.MinimumEfficiency");
  Param.GetScalar(SubgridSizeAutoAdjust, "PhysicsParameters.Optimization.SubgridSizeAutoAdjust");
  Param.GetScalar(OptimalSubgridsPerProcessor, "PhysicsParameters.Optimization.OptimalSubgridsPerProcessor");
  Param.GetScalar(MinimumSubgridEdge, "PhysicsParameters.Optimization.MinimumSubgridEdge");
  Param.GetScalar(MaximumSubgridSize, "PhysicsParameters.Optimization.MaximumSubgridSize");
  Param.GetScalar(NumberOfBufferZones, "SimulationControl.AMR.NumberOfBufferZones");
  Param.GetScalar(FastSiblingLocatorEntireDomain, "PhysicsParameters.Optimization.FastSiblingLocatorEntireDomain");
  Param.GetScalar(MustRefineRegionMinRefinementLevel, "SimulationControl.AMR.MustRefineRegionMinRefinementLevel");
  Param.GetScalar(MetallicityRefinementMinLevel, "SimulationControl.AMR.MetallicityRefinementMinLevel");  
  Param.GetScalar(MetallicityRefinementMinMetallicity, "SimulationControl.AMR.MetallicityRefinementMinMetallicity"); 
  Param.GetScalar(MetallicityRefinementMinDensity, "SimulationControl.AMR.MetallicityRefinementMinDensity"); 
  Param.GetArray(DomainLeftEdge, "SimulationControl.Domain.DomainLeftEdge");
  Param.GetArray(DomainRightEdge, "SimulationControl.Domain.DomainRightEdge");
  Param.GetArray(GridVelocity, "PhysicsParameters.GridVelocity");
  
  Param.GetScalar(RefineRegionAutoAdjust, "SimulationControl.Domain.RefineRegionAutoAdjust"); 
  Param.GetArray(RefineRegionLeftEdge, "SimulationControl.Domain.RefineRegionLeftEdge");
  Param.GetArray(RefineRegionRightEdge, "SimulationControl.Domain.RefineRegionRightEdge");
  Param.GetArray(MustRefineRegionLeftEdge, "SimulationControl.Domain.MustRefineRegionLeftEdge");
  Param.GetArray(MustRefineRegionRightEdge, "SimulationControl.Domain.MustRefineRegionRightEdge");
  
  
  /* Read evolving RefineRegion */
  
  Param.GetScalar(RefineRegionTimeType, "SimulationControl.AMR.RefineRegionTimeType");
  Param.GetScalar(RefineRegionFile, "SimulationControl.AMR.RefineRegionFile");
  
  int NumberOfFields = Param.Size("InternalParameters.Fields");
  if (NumberOfFields > MAX_NUMBER_OF_BARYON_FIELDS) {
    ENZO_VFAIL("You've exceeded the maximum number of BaryonFields (%d)!\n",MAX_NUMBER_OF_BARYON_FIELDS)
      }
  
  char FieldNames[MAX_LINE_LENGTH][NumberOfFields];
  Param.GetArray(StaticRefineRegionNames, "InternalParameters.Fields");
  
  for (i = 0; i < NumberOfFields; i++) {
    Param.GetScalar(DataLabel[i], "InternalParameters.%S.Name", FieldNames[i]);
    Param.GetScalar(DataUnits[i], "InternalParameters.%s.cgsConversionFactor", FieldNames[i]);
  }
  
  Param.GetScalar(UniformGravity, "PhysicsParameters.Gravity.UniformGravity");
  Param.GetScalar(UniformGravityDirection, "PhysicsParameters.Gravity.UniformGravityDirection");
  Param.GetScalar(UniformGravityConstant, "PhysicsParameters.Gravity.UniformGravityConstant");
  
  Param.GetScalar(PointSourceGravity, "PhysicsParameters.Gravity.PointSourceGravity");
  Param.GetArray(PointSourceGravityPosition, "PhysicsParameters.Gravity.PointSourceGravityPosition");
  Param.GetScalar(PointSourceGravityConstant, "PhysicsParameters.Gravity.PointSourceGravityConstant");
  Param.GetScalar(PointSourceGravityCoreRadius, "PhysicsParameters.Gravity.PointSourceGravityCoreRadius");
  
  Param.GetScalar(ExternalGravity, "PhysicsParameters.Gravity.ExternalGravity");
  
  Param.GetScalar(SelfGravity, "PhysicsParameters.Gravity.SelfGravity");
  Param.GetScalar(SelfGravityGasOff, "PhysicsParameters.Gravity.SelfGravityGasOff");
  Param.GetScalar(AccretionKernel, "PhysicsParameters.AccretionKernel");
  Param.GetScalar(GravitationalConstant, "PhysicsParameters.Gravity.GravitationalConstant");
  
  Param.GetScalar(S2ParticleSize, "PhysicsParameters.OtherParticles.S2ParticleSize");
  Param.GetScalar(GravityResolution, "PhysicsParameters.Gravity.GravityResolution");
  Param.GetScalar(ComputePotential, "PhysicsParameters.Gravity.ComputePotential");
  Param.GetScalar(PotentialIterations, "PhysicsParameters.Gravity.PotentialIterations");
  Param.GetScalar(WritePotential, "OutputControlParameters.SupplementalFields.WritePotential"); 
  Param.GetScalar(BaryonSelfGravityApproximation, "PhysicsParameters.Gravity.BaryonSelfGravityApproximation"); 
  
  Param.GetScalar(GreensFunctionMaxNumber, "PhysicsParameters.Gravity.GreensFunctionMaxNumber");
  Param.GetScalar(GreensFunctionMaxSize, "PhysicsParameters.Gravity.GreensFunctionMaxSize");
  
  Param.GetScalar(DualEnergyFormalism, "PhysicsParameters.Hydro.DualEnergyFormalism");
  Param.GetScalar(DualEnergyFormalismEta1, "PhysicsParameters.Hydro.DualEnergyFormalismEta1");
  Param.GetScalar(DualEnergyFormalismEta2, "PhysicsParameters.Hydro.DualEnergyFormalismEta2");
  Param.GetScalar(ParticleCourantSafetyNumber, "PhysicsParameters.Hydro.ParticleCourantSafetyNumber");
  Param.GetScalar(RootGridCourantSafetyNumber, "PhysicsParameters.Hydro.RootGridCourantSafetyNumber");
  Param.GetScalar(RandomForcing, "PhysicsParameters.Miscellaneous.RandomForcing");
  Param.GetScalar(RandomForcingEdot, "PhysicsParameters.Miscellaneous.RandomForcingEdot");
  Param.GetScalar(RandomForcingMachNumber, "PhysicsParameters.Miscellaneous.RandomForcingMachNumber");
  Param.GetScalar(RadiativeCooling, "PhysicsParameters.AtomicPhysics.RadiativeCooling");
  Param.GetScalar(RadiativeCoolingModel, "PhysicsParameters.AtomicPhysics.RadiativeCoolingModel");
  Param.GetScalar(GadgetEquilibriumCooling, "PhysicsParameters.AtomicPhysics.GadgetEquilibriumCooling");
  Param.GetScalar(MultiSpecies, "PhysicsParameters.AtomicPhysics.MultiSpecies");
  Param.GetScalar(PrimordialChemistrySolver, "PhysicsParameters.AtomicPhysics.PrimordialChemistrySolver");
  Param.GetScalar(CIECooling, "PhysicsParameters.AtomicPhysics.CIECooling");
  Param.GetScalar(H2OpticalDepthApproximation, "PhysicsParameters.AtomicPhysics.H2OpticalDepthApproximation");
  Param.GetScalar(ThreeBodyRate, "PhysicsParameters.AtomicPhysics.ThreeBodyRate");
  
  
  Param.GetScalar(CloudyCoolingData.CloudyCoolingGridFile, "PhysicsParameters.AtomicPhysics.CloudyCooling");
  Param.GetScalar(CloudyCoolingData.IncludeCloudyHeating, "PhysicsParameters.AtomicPhysics.CloudyCooling.IncludeCloudyHeating");
  Param.GetScalar(CloudyCoolingData.IncludeCloudyMMW, "PhysicsParameters.AtomicPhysics.CloudyCooling.IncludeCloudyMMW");
  Param.GetScalar(CloudyCoolingData.CMBTemperatureFloor, "PhysicsParameters.AtomicPhysics.CloudyCooling.CMBTemperatureFloor");
  Param.GetScalar(CloudyCoolingData.CloudyMetallicityNormalization, "PhysicsParameters.AtomicPhysics.CloudyCooling.CloudyMetallicityNormalization");
  Param.GetScalar(CloudyCoolingData.CloudyElectronFractionFactor, "PhysicsParameters.AtomicPhysics.CloudyCooling.CloudyElectronFractionFactor");
  Param.GetScalar(MetalCooling, "PhysicsParameters.AtomicPhysics.MetalCooling");
  
  Param.GetScalar(MetalCoolingTable, "PhysicsParameters.AtomicPhysics.MetalCoolingTable");
  
  Param.GetScalar(CRModel, "PhysicsParameters.Miscellaneous.CRModel");
  Param.GetScalar(ShockMethod, "PhysicsParameters.Miscellaneous.ShockMethod");
  Param.GetScalar(ShockTemperatureFloor, "PhysicsParameters.Miscellaneous.ShockTemperatureFloor");
  Param.GetScalar(StorePreShockFields, "PhysicsParameters.Miscellaneous.StorePreShockFields");
  
  Param.GetScalar(RadiationFieldType, "PhysicsParameters.RadiationField.RadiationFieldType");
  Param.GetScalar(TabulatedLWBackground, "PhysicsParameters.RadiationField.TabulatedLWBackground");
  Param.GetScalar(AdjustUVBackground, "PhysicsParameters.RadiationField.AdjustUVBackground");
  Param.GetScalar(SetUVBAmplitude, "PhysicsParameters.RadiationField.SetUVBAmplitude");
  Param.GetScalar(SetHeIIHeatingScale, "PhysicsParameters.RadiationField.SetHeIIHeatingScale");
  Param.GetScalar(RadiationFieldLevelRecompute, "PhysicsParameters.RadiationField.RadiationFieldLevelRecompute");
  Param.GetScalar(RadiationData.RadiationShield, "PhysicsParameters.RadiationField.RadiationShield");
  Param.GetScalar(CoolData.f3, "PhysicsParameters.RadiationField.RadiationSpectrumNormalization");
  Param.GetScalar(CoolData.alpha0, "PhysicsParameters.RadiationField.RadiationSpectrumSlope");
  Param.GetScalar(CoolData.f0to3, "PhysicsParameters.AtomicPhysics.CoolDataf0to3");
  Param.GetScalar(CoolData.RadiationRedshiftOn, "PhysicsParameters.RadiationField.RadiationRedshiftOn");
  Param.GetScalar(CoolData.RadiationRedshiftOff, "PhysicsParameters.RadiationField.RadiationRedshiftOff");
  Param.GetScalar(CoolData.RadiationRedshiftFullOn, "PhysicsParameters.RadiationField.RadiationRedshiftFullOn");
  Param.GetScalar(CoolData.RadiationRedshiftDropOff, "PhysicsParameters.RadiationField.RadiationRedshiftDropOff");
  Param.GetScalar(CoolData.HydrogenFractionByMass, "PhysicsParameters.AtomicPhysics.HydrogenFractionByMass");
  Param.GetScalar(CoolData.DeuteriumToHydrogenRatio, "PhysicsParameters.AtomicPhysics.DeuteriumToHydrogenRatio");
  Param.GetScalar(CoolData.NumberOfTemperatureBins, "PhysicsParameters.AtomicPhysics.NumberOfTemperatureBins");
  Param.GetScalar(CoolData.ih2co, "PhysicsParameters.AtomicPhysics.CoolDataIh2co");
  Param.GetScalar(CoolData.ipihdt, "PhysicsParameters.AtomicPhysics.CoolDataIpiht");
  Param.GetScalar(CoolData.TemperatureStart, "PhysicsParameters.AtomicPhysics.TemperatureStart");
  Param.GetScalar(CoolData.TemperatureEnd, "PhysicsParameters.AtomicPhysics.TemperatureEnd");
  Param.GetScalar(CoolData.comp_xray, "PhysicsParameters.AtomicPhysics.CoolDataCompXray");
  Param.GetScalar(CoolData.temp_xray, "PhysicsParameters.AtomicPhysics.CoolDataTempXray");
  Param.GetScalar(RateData.CaseBRecombination, "PhysicsParameters.AtomicPhysics.RateDataCaseBRecombination");
  Param.GetScalar(PhotoelectricHeating, "PhysicsParameters.AtomicPhysics.PhotoelectricHeating");
  
  Param.GetScalar(OutputCoolingTime, "OutputControlParameters.SupplementalFields.OutputCoolingTime"); 
  Param.GetScalar(OutputTemperature, "OutputControlParameters.SupplementalFields.OutputTemperature"); 
  Param.GetScalar(OutputSmoothedDarkMatter, "OutputControlParameters.SupplementalFields.OutputSmoothedDarkMatter"); 
  Param.GetScalar(SmoothedDarkMatterNeighbors, "OutputControlParameters.SupplementalFields.SmoothedDarkMatterNeighbors"); 
  Param.GetScalar(OutputGriddedStarParticle, "OutputControlParameters.SupplementalFields.OutputGriddedStarParticle"); 
  
  Param.GetScalar(ZEUSQuadraticArtificialViscosity, "PhysicsParameters.Hydro.ZEUSQuadraticArtificialViscosity");
  Param.GetScalar(ZEUSLinearArtificialViscosity, "PhysicsParameters.Hydro.ZEUSLinearArtificialViscosity");
  
  Param.GetScalar(UseMinimumPressureSupport, "PhysicsParameters.Hydro.UseMinimumPressureSupport"); 
  Param.GetScalar(MinimumPressureSupportParameter, "PhysicsParameters.Hydro.MinimumPressureSupportParameter");
  Param.GetScalar(RefineByJeansLengthSafetyFactor, "SimulationControl.AMR.RefineByJeansLengthSafetyFactor");
  Param.GetScalar(RefineByJeansLengthUnits, "SimulationControl.AMR.RefineByJeansLengthUnits"); 
  Param.GetScalar(JeansRefinementColdTemperature, "SimulationControl.AMR.JeansRefinementColdTemperature");
  Param.GetScalar(RefineByResistiveLengthSafetyFactor, "SimulationControl.AMR.RefineByResistiveLengthSafetyFactor");
  Param.GetScalar(MustRefineParticlesRefineToLevel, "SimulationControl.AMR.MustRefineParticlesRefineToLevel");
  Param.GetScalar(MustRefineParticlesRefineToLevelAutoAdjust, "SimulationControl.AMR.MustRefineParticlesRefineToLevelAutoAdjust"); 
  Param.GetScalar(MustRefineParticlesMinimumMass, "SimulationControl.AMR.MustRefineParticlesMinimumMass");
  Param.GetScalar(ParticleTypeInFile, "OutputControlParameters.ParticleTypeInFile"); 
  
  int NumberOfStaticRefineRegions = Param.Size("StaticRefineRegion.Regions");
  if (NumberOfStaticRefineRegions > MAX_STATIC_REGIONS-1) {
    ENZO_VFAIL("You've exceeded the maximum number of StaticRefineRegions (%d)!\n",MAX_STATIC_REGIONS)
      }
  
  char StaticRefineRegionNames[MAX_LINE_LENGTH][NumberOfStaticRefineRegions];
  Param.GetArray(StaticRefineRegionNames, "StaticRefineRegion.Regions");
  
  for (i = 0; i < NumberOfStaticRefineRegions; i++) {
    Param.GetScalar(StaticRefineRegionLevel[i], "StaticRefineRegion.%s.Level",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionLeftEdge[i],"StaticRefineRegion.%s.LeftEdge",StaticRefineRegionNames[i]);
    Param.GetArray(StaticRefineRegionRightEdge[i],"StaticRefineRegion.%s.RightEdge",StaticRefineRegionNames[i]);
  }
  
  Param.GetScalar(ParallelRootGridIO, "InternalParameters.OutputLabeling.ParallelRootGridIO");
  Param.GetScalar(ParallelParticleIO, "InternalParameters.OutputLabeling.ParallelParticleIO");
  
  Param.GetScalar(Unigrid, "InternalParameters.OutputLabeling.Unigrid");
  Param.GetScalar(UnigridTranspose, "InternalParameters.OutputLabeling.UnigridTranspose");
  Param.GetScalar(NumberOfRootGridTilesPerDimensionPerProcessor, "InternalParameters.OutputLabeling.NumberOfRootGridTilesPerDimensionPerProcessor");
  
  Param.GetScalar(PartitionNestedGrids, "Initialization.PartitionNestedGrids");
  Param.GetScalar(StaticPartitionNestedGrids, "Initialization.StaticPartitionNestedGrids");
  
  Param.GetScalar(ExtractFieldsOnly, "OutputControlParameters.ExtractFieldsOnly");
  
  Param.GetScalar(debug1, "InternalParameters.Debug1");
  
  Param.GetScalar(debug2, "InternalParameters.Debug2");
  
  Param.GetScalar(MemoryLimit, "SimulationControl.Optimization.MemoryLimit");
  
  
#ifdef OOC_BOUNDARY
  
  Param.GetScalar(ExternalBoundaryIO, "OutputControlParameters.SupplementalFields.ExternalBoundaryIO");
  
  Param.GetScalar(ExternalBoundaryTypeIO, "OutputControlParameters.SupplementalFields.ExternalBoundaryTypeIO");
  
  Param.GetScalar(ExternalBoundaryValueIO, "OutputControlParameters.SupplementalFields.ExternalBoundaryValueIO");
  
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
  Param.GetScalar(ComovingCoordinates, "PhysicsParameters.Cosmology.ComovingCoordinates");
  Param.GetScalar(StarParticleCreation, "PhysicsParameters.OtherParticles.StarParticleCreation");
  Param.GetScalar(BigStarFormation, "PhysicsParameters.OtherParticles.BigStar.BigStarFormation");
  Param.GetScalar(BigStarFormationDone, "PhysicsParameters.OtherParticles.BigStar.BigStarFormationDone");
  Param.GetScalar(BigStarSeparation, "PhysicsParameters.OtherParticles.BigStar.BigStarSeparation");
  Param.GetScalar(SimpleQ, "SimulationControl.RadiationTransfer.SimpleQ");
  Param.GetScalar(SimpleRampTime, "SimulationControl.RadiationTransfer.SimpleRampTime");
  Param.GetScalar(StarParticleFeedback, "PhysicsParameters.OtherParticles.StarParticleFeedback");
  Param.GetScalar(NumberOfParticleAttributes, "PhysicsParameters.OtherParticles.NumberOfParticleAttributes");
  Param.GetScalar(AddParticleAttributes, "PhysicsParameters.OtherParticles.AddParticleAttributes");
  
  
  /* read data which defines the boundary conditions */
  Param.GetArray(MetaData.LeftFaceBoundaryCondition, "PhysicsParameters.LeftFaceBoundaryCondition");
  Param.GetArray(MetaData.RightFaceBoundaryCondition, "PhysicsParameters.RightFaceBoundaryCondition");
  
  Param.GetScalar(MetaData.BoundaryConditionName, "InternalParameters.BoundaryConditionName");
  Param.GetScalar(MetaData.MetaDataIdentifier, "InternalParameters.Provenance.MetaDataIdentifier");
  Param.GetScalar(MetaData.SimulationUUID, "InternalParameters.Provenance.MetaDataSimulationUUID");
  Param.GetScalar(MetaData.RestartDatasetUUID, "InternalParameters.Provenance.MetaDataDatasetUUID");
  Param.GetScalar(MetaData.InitialConditionsUUID, "InternalParameters.Provenance.MetaDataInitialConditionsUUID");
  
  /* Check version number. */
  
  Param.GetScalar(TempFloat, "InternalParameters.Provenance.VersionNumber");
  if (fabs(TempFloat - VERSION) >= 1.0e-3 &&
      MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Warning: Incorrect version number.\n");
  
  
  /* Read star particle parameters. */
  
  Param.GetScalar(StarMakerOverDensityThreshold, "PhysicsParameters.OtherParticles.StarMaker.StarMakerOverDensityThreshold");
  Param.GetScalar(StarMakerSHDensityThreshold, "PhysicsParameters.OtherParticles.StarMaker.StarMakerSHDensityThreshold");
  Param.GetScalar(StarMakerMassEfficiency, "PhysicsParameters.OtherParticles.StarMaker.StarMakerMassEfficiency");
  Param.GetScalar(StarMakerMinimumMass, "PhysicsParameters.OtherParticles.StarMaker.StarMakerMinimumMass");
  Param.GetScalar(StarMakerMinimumDynamicalTime, "PhysicsParameters.OtherParticles.StarMaker.StarMakerMinimumDynamicalTime");
  Param.GetScalar(StarMassEjectionFraction, "PhysicsParameters.OtherParticles.StarMaker.StarMassEjectionFraction");
  Param.GetScalar(StarMetalYield, "PhysicsParameters.OtherParticles.StarMaker.StarMetalYield");
  Param.GetScalar(StarEnergyToThermalFeedback, "PhysicsParameters.OtherParticles.StarEnergy.StarEnergyToThermalFeedback");
  Param.GetScalar(StarEnergyToStellarUV, "PhysicsParameters.OtherParticles.StarEnergy.StarEnergy.ToStellarUV");
  Param.GetScalar(StarEnergyToQuasarUV, "PhysicsParameters.OtherParticles.StarEnergy.StarEnergyToQuasarUV");
  
  Param.GetScalar(StarFeedbackDistRadius, "PhysicsParameters.OtherParticles.StarFeedbackDistRadius");
  Param.GetScalar(StarFeedbackDistCellStep, "PhysicsParameters.OtherParticles.StarFeedbackDistCellStep");
  Param.GetScalar(StarClusterUseMetalField, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterUseMetalField");
  Param.GetScalar(StarClusterMinDynamicalTime, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterMinDynamicalTime");
  Param.GetScalar(StarClusterHeliumIonization, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterHeliumIonization");
  Param.GetScalar(StarClusterUnresolvedModel, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterUnresolvedModel");
  Param.GetScalar(StarClusterIonizingLuminosity, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterIonizingLuminosity");
  Param.GetScalar(StarClusterSNEnergy, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterSNEnergy");
  Param.GetScalar(StarClusterSNRadius, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterSNRadius");
  Param.GetScalar(StarClusterFormEfficiency, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterFormEfficiency");
  Param.GetScalar(StarClusterMinimumMass, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterMinimumMass");
  Param.GetScalar(StarClusterCombineRadius, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterCombineRadius");
  
  Param.GetArray(StarClusterRegionLeftEdge, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterRegionLeftEdge");
  Param.GetArray(StarClusterRegionRightEdge, "PhysicsParameters.OtherParticles.StarClusterParticle.StarClusterRegionRightEdge");
  
  Param.GetScalar(PopIIIStarMass, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIStarMass");                           
  Param.GetScalar(PopIIIInitialMassFunctionSeed, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIInitialMassFunctionSeed");                 
  Param.GetScalar(PopIIIInitialMassFunctionCalls, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIInitialMassFunctionCalls"); 
  Param.GetScalar(PopIIILowerMassCutoff, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIILowerMassCutoff");
  Param.GetScalar(PopIIIUpperMassCutoff, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIUpperMassCutoff");
  
  //   ret += sscanf(line,"PopIIIMassRange %"FSYM" %"FSYM, &PopIIILowerMassCutoff, &PopIIIUpperMassCutoff);
  
  Param.GetScalar(PopIIIInitialMassFunctionSlope, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIInitialMassFunctionSlope");
  Param.GetScalar(PopIIIBlackHoles, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIBlackHoles");
  Param.GetScalar(PopIIIBHLuminosityEfficiency, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIBHLuminosityEfficiency");
  Param.GetScalar(PopIIIOverDensityThreshold, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIOverDensityThreshold");
  Param.GetScalar(PopIIIH2CriticalFraction, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIH2CriticalFraction");
  Param.GetScalar(PopIIIMetalCriticalFraction, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIMetalCriticalFraction");
  Param.GetScalar(PopIIISupernovaRadius, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIISupernovaRadius");
  Param.GetScalar(PopIIISupernovaUseColour, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIISupernovaUseColour");
  Param.GetScalar(PopIIISupernovaMustRefine, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIISupernovaMustRefine");
  Param.GetScalar(PopIIISupernovaMustRefineResolution, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIISupernovaMustRefineResolution");
  Param.GetScalar(PopIIIHeliumIonization, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIHeliumIonization");
  Param.GetScalar(PopIIIColorDensityThreshold, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIColorDensityThreshold");
  Param.GetScalar(PopIIIColorMass, "PhysicsParameters.OtherParticles.PopIIIStarParticle.PopIIIColorMass");
  
  Param.GetScalar(MBHAccretion, "PhysicsParameters.OtherParticles.MBHParticle.MBHAccretion");
  Param.GetScalar(MBHAccretionRadius, "PhysicsParameters.OtherParticles.MBHParticle.MBHAccretionRadius");
  Param.GetScalar(MBHAccretingMassRatio, "PhysicsParameters.OtherParticles.MBHParticle.MBHAccretingMassRatio");
  Param.GetScalar(MBHAccretionFixedTemperature, "PhysicsParameters.OtherParticles.MBHParticle.MBHAccretionFixedTemperature");
  Param.GetScalar(MBHAccretionFixedRate, "PhysicsParameters.OtherParticles.MBHParticle.MBHAccretionFixedRate");
  Param.GetScalar(MBHTurnOffStarFormation, "PhysicsParameters.OtherParticles.MBHParticle.MBHTurnOffStarFormation");
  Param.GetScalar(MBHCombineRadius, "PhysicsParameters.OtherParticles.MBHParticle.MBHCombineRadius");
  Param.GetScalar(MBHMinDynamicalTime, "PhysicsParameters.OtherParticles.MBHParticle.MBHMinDynamicalTime");
  Param.GetScalar(MBHMinimumMass, "PhysicsParameters.OtherParticles.MBHParticle.MBHMinimumMass");
  
  Param.GetScalar(MBHFeedback, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedback");
  Param.GetScalar(MBHFeedbackRadiativeEfficiency, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackRadiativeEfficiency");
  Param.GetScalar(MBHFeedbackEnergyCoupling, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackEnergyCoupling");
  Param.GetScalar(MBHFeedbackMassEjectionFraction, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackMassEjectionFraction");
  Param.GetScalar(MBHFeedbackMetalYield, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackMetalYield");
  Param.GetScalar(MBHFeedbackThermalRadius, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackThermalRadius");
  Param.GetScalar(MBHFeedbackJetsThresholdMass, "PhysicsParameters.OtherParticles.MBHParticle.MBHFeedbackJetsThresholdMass");
  Param.GetScalar(MBHParticleIO, "PhysicsParameters.OtherParticles.MBHParticle.MBHParticleIO");
  
  Param.GetScalar(MBHParticleIOFilename, "Initialization.MBHParticleIOFilename");
  Param.GetScalar(MBHInsertLocationFilename, "Initialization.MBHInsertLocationFilename");
  
  
  /* Read Movie Dump parameters */
  
  Param.GetScalar(MovieSkipTimestep, "OutputControlParameters.MovieDump.MovieSkipTimestep");
  Param.GetScalar(Movie3DVolumes, "OutputControlParameters.MovieDump.Movie3DVolumes");
  Param.GetScalar(MovieVertexCentered, "OutputControlParameters.MovieDump.MovieVertexCentered");
  Param.GetScalar(NewMovieParticleOn, "OutputControlParameters.MovieDump.NewMovieParticleOn");
  Param.GetArray(MovieDataField, "OutputControlParameters.MovieDump.MovieDataField");
  Param.GetScalar(NewMovieDumpNumber, "OutputControlParameters.MovieDump.NewMovieDumpNumber");
  Param.GetScalar(NewMovieName, "OutputControlParameters.MovieDump.NewMovieName");
  Param.GetScalar(MetaData.MovieTimestepCounter, "OutputControlParameters.MovieDump.MovieTimestepCounter");
  
  Param.GetScalar(MultiMetals, "PhysicsParameters.AtomicPhysics.MultiMetals");
  Param.GetScalar(IsotropicConduction, "PhysicsParameters.Conduction.IsotropicConduction");
  Param.GetScalar(AnisotropicConduction, "PhysicsParameters.Conduction.AnisotropicConduction");
  Param.GetScalar(IsotropicConductionSpitzerFraction, "PhysicsParameters.Conduction.IsotropicConductionSpitzerFraction");
  Param.GetScalar(AnisotropicConductionSpitzerFraction, "PhysicsParameters.Conduction.AnisotropicConductionSpitzerFraction");
  Param.GetScalar(ConductionCourantSafetyNumber, "PhysicsParameters.Conduction.ConductionCourantSafetyNumber");
  
  Param.GetScalar(RadiativeTransfer, "PhysicsParameters.RadiationTransfer.RadiativeTransfer");
  Param.GetScalar(RadiationXRaySecondaryIon, "PhysicsParameters.RadiationTransfer.RadiationXRaySecondaryIon");
  Param.GetScalar(RadiationXRayComptonHeating, "PhysicsParameters.RadiationTransfer.RadiationXRayComptonHeating");
  
  
  /* Shearing Box Boundary parameters */
  
  Param.GetScalar(AngularVelocity, "PhysicsParameters.AngularVelocity");
  Param.GetScalar(VelocityGradient, "PhysicsParameters.VelocityGradient");
  Param.GetScalar(ShearingVelocityDirection, "PhysicsParameters.ShearingVelocityDirection");
  Param.GetScalar(ShearingBoxProblemType, "PhysicsParameters.ShearingBoxProblemType");
  
  /* Embedded Python */
  Param.GetScalar(PythonTopGridSkip, "AnalysisParameters.Python.PythonTopGridSkip");
  Param.GetScalar(PythonSubcycleSkip, "AnalysisParameters.Python.PythonSubcycleSkip");

#ifdef USE_PYTHON
  Param.GetScalar(NumberOfPythonCalls, "AnalysisParameters.Python.NumberOfPythonCalls");
  Param.GetScalar(NumberOfPythonTopGridCalls, "AnalysisParameters.Python.NumberOfPythonTopGridCalls");
  Param.GetScalar(NumberOfPythonSubcycleCalls, "AnalysisParameters.Python.NumberOfPythonSubcycleCalls");
#endif

    /* Inline halo finder */

  Param.GetScalar(InlineHaloFinder, "AnalysisParameters.HaloFinder.InlineHaloFinder");
  Param.GetScalar(HaloFinderSubfind, "AnalysisParameters.HaloFinder.HaloFinderSubfind");
  Param.GetScalar(HaloFinderOutputParticleList, "AnalysisParameters.HaloFinder.HaloFinderOutputParticleList");
  Param.GetScalar(HaloFinderRunAfterOutput, "AnalysisParameters.HaloFinder.HaloFinderRunAfterOutput");
  Param.GetScalar(HaloFinderLinkingLength, "AnalysisParameters.HaloFinder.HaloFinderLinkingLength");
  Param.GetScalar(HaloFinderMinimumSize, "AnalysisParameters.HaloFinder.HaloFinderMinimumSize");
  Param.GetScalar(HaloFinderCycleSkip, "AnalysisParameters.HaloFinder.HaloFinderCycleSkip");
  Param.GetScalar(HaloFinderTimestep, "AnalysisParameters.HaloFinder.HaloFinderTimestep");
  Param.GetScalar(HaloFinderLastTime, "AnalysisParameters.HaloFinder.HaloFinderLastTime");
  
  /* This Block for Stanford Hydro */
  
  Param.GetScalar(UseHydro, "PhysicsParameters.Hydro.UseHydro");
  
  
  /* Sink particles (for present day star formation) & winds */
  Param.GetScalar(SinkMergeDistance, "PhysicsParameters.OtherParticles.SinkParticle.SinkMergeDistance");
  Param.GetScalar(SinkMergeMass, "PhysicsParameters.OtherParticles.SinkParticle.SinkMergeMass");
  Param.GetScalar(StellarWindFeedback, "PhysicsParameters.OtherParticles.StellarWindFeedback");
  Param.GetScalar(StellarWindTurnOnMass, "PhysicsParameters.OtherParticles.StellarWindTurnOnMass");
  Param.GetScalar(MSStellarWindTurnOnMass, "PhysicsParameters.OtherParticles.MSStellarWindTurnOnMass");
  
  Param.GetScalar(VelAnyl, "OutputControlParameters.SupplementalFields.VelAnyl");
  Param.GetScalar(BAnyl, "OutputControlParameters.SupplementalFields.BAnyl");
  
  
  
  /* Read MHD Paramters */
  Param.GetScalar(UseDivergenceCleaning, "PhysicsParameters.MHD.UseDivergenceCleaning");
  Param.GetScalar(DivergenceCleaningBoundaryBuffer, "PhysicsParameters.MHD.DivergenceCleaningBoundaryBuffer");
  Param.GetScalar(DivergenceCleaningThreshold, "PhysicsParameters.MHD.DivergenceCleaningThreshold");
  Param.GetScalar(PoissonApproximationThreshold, "PhysicsParameters.PoissonApproximationThreshold");
  Param.GetScalar(PoissonBoundaryType, "PhysicsParameters.PoissonBoundaryType");
  
  
  Param.GetScalar(AngularVelocity, "PhysicsParameters.AngularVelocity");
  Param.GetScalar(VelocityGradient, "PhysicsParameters.VelocityGradient");
  Param.GetScalar(UseDrivingField, "PhysicsParameters.Miscellaneous.UseDrivingField");
  Param.GetScalar(DrivingEfficiency, "PhysicsParameters.Miscellaneous.DrivingEfficiency");
  
  Param.GetScalar(StringKick, "Initialization.StringKick");
  Param.GetScalar(StringKickDimension, "Initialization.StringKickDimension");
  Param.GetScalar(UsePhysicalUnit, "PhysicsParameters.UsePhysicalUnit");
  Param.GetScalar(Theta_Limiter, "PhysicsParameters.Hydro.Theta_Limiter");
  Param.GetScalar(RKOrder, "PhysicsParameters.Hydro.RKOrder");
  Param.GetScalar(UseFloor, "PhysicsParameters.Hydro.UseFloor");
  Param.GetScalar(UseViscosity, "PhysicsParameters.Hydro.UseViscosity");
  Param.GetScalar(ViscosityCoefficient, "PhysicsParameters.Hydro.ViscosityCoefficient");
  
  Param.GetScalar(UseAmbipolarDiffusion, "PhysicsParameters.Hydro.UseAmbipolarDiffusion");
  Param.GetScalar(UseResistivity, "PhysicsParameters.Hydro.UseResistivity");
  Param.GetScalar(SmallRho, "PhysicsParameters.SmallRho");
  Param.GetScalar(SmallP, "PhysicsParameters.SmallP");
  Param.GetScalar(SmallT, "PhysicsParameters.SmallT");
  Param.GetScalar(MaximumAlvenSpeed, "PhysicsParameters.MHD.MaximumAlvenSpeed");
  Param.GetScalar(Coordinate, "PhysicsParameters.Coordinate");
  Param.GetScalar(RiemannSolver, "PhysicsParameters.Hydro.RiemannSolver");
  Param.GetScalar(RiemannSolverFallback, "PhysicsParameters.Hydro.RiemannSolverFallback");
  Param.GetScalar(ConservativeReconstruction, "PhysicsParameters.Hydro.ConservativeReconstruction");
  Param.GetScalar(PositiveReconstruction, "PhysicsParameters.Hydro.PositiveReconstruction");
  Param.GetScalar(ReconstructionMethod, "PhysicsParameters.Hydro.ReconstructionMethod");
  
  Param.GetScalar(EOSType, "PhysicsParameters.Hydro.EOSType");
  Param.GetScalar(EOSSoundSpeed, "PhysicsParameters.Hydro.EOSSoundSpeed");
  Param.GetScalar(EOSCriticalDensity, "PhysicsParameters.Hydro.EOSCriticalDensity");
  Param.GetScalar(EOSGamma, "PhysicsParameters.Hydro.EOSGamma");
  Param.GetScalar(UseConstantAcceleration, "PhysicsParameters.Gravity.UseConstantAcceleration");
  
  Param.GetScalar(IsothermalSoundSpeed, "PhysicsParameters.Hydro.IsothermalSoundSpeed");
  
  Param.GetArray(ConstantAcceleration, "PhysicsParameters.Gravity.ConstantAcceleration");
  Param.GetScalar(Mu, "PhysicsParameters.Hydro.Mu");
  Param.GetScalar(DivBDampingLength, "PhysicsParameters.MHD.DivBDampingLength");
  Param.GetScalar(CoolingCutOffDensity1, "PhysicsParameters.AtomicPhysics.CoolingCutOffDensity1");
  Param.GetScalar(CoolingCutOffDensity2, "PhysicsParameters.AtomicPhysics.CoolingCutOffDensity2");
  Param.GetScalar(CoolingCutOffTemperature, "PhysicsParameters.AtomicPhysics.CoolingCutOffTemperature");
  Param.GetScalar(CoolingPowerCutOffDensity1, "PhysicsParameters.AtomicPhysics.CoolingPowerCutOffDensity1");
  Param.GetScalar(CoolingPowerCutOffDensity2, "PhysicsParameters.AtomicPhysics.CoolingPowerCutOffDensity2");
  Param.GetScalar(UseH2OnDust, "PhysicsParameters.AtomicPhysics.UseH2OnDust");
  Param.GetScalar(seCUDA, "SimulationControl.UseCUDA");
  
  Param.GetScalar(MoveParticlesBetweenSiblings, "SimulationControl.Optimization.MoveParticlesBetweenSiblings");
  Param.GetScalar(ParticleSplitterIterations, "SimulationControl.Optimization.ParticleSplitterIterations");
  Param.GetScalar(ParticleSplitterChildrenParticleSeparation, "SimulationControl.Optimization.ParticleSplitterChildrenParticleSeparation");
  Param.GetScalar(ResetMagneticField, "PhysicsParameters.MHD.ResetMagneticField");
  Param.GetArray(ResetMagneticFieldAmplitude, "PhysicsParameters.MHD.ResetMagneticFieldAmplitude");

  
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
