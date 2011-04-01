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
// #include <string.h>
#include <string>
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

int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt)
{
#ifndef CONFIG_USE_LIBCONFIG
  /* declarations */

  
  char line[MAX_LINE_LENGTH];
  int i, dim, ret, int_dummy;
  float TempFloat;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  int comment_count = 0;
 
  /* read until out of lines */

  rewind(fptr);
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {

    ret = 0;
 
    /* read MetaData parameters */
 
    MetaData.CycleNumber = Param.GetScalar <int> ("InternalParameters.InitialCycleNumber");
    MetaData.Time        = Param.GetScalar <FLOAT> ("InternalParameters.InitialTime");
    MetaData.CPUTime     = Param.GetScalar <float> ("InternalParameters.InitialCPUTime");
    (*Initialdt)         = Param.GetScalar <float> ("InternalParameters.Initialdt");
 
    CheckpointRestart    = Param.GetScalar <int> ("CheckpointRestart"); // should be bool
    MetaData.StopTime    = Param.GetScalar <FLOAT> ("SimulationControl.StopTime");
    MetaData.StopCycle   = Param.GetScalar <int> ("SimulationControl.StopCycle");
    MetaData.StopSteps   = Param.GetScalar <int> ("SimulationControl.StopSteps");
    MetaData.StopCPUTime = Param.GetScalar <float> ("StopCPUTime");
    MetaData.ResubmitOn  = Param.GetScalar <int> ("SimulationControl.ResubmitOn"); // should be bool
    /* check syntex of sscanf to MetaData */
    if (sscanf(line, "ResubmitCommand = %s", dummy) == 1) MetaData.ResubmitCommand = dummy;

    MetaData.MaximumTopGridTimeStep = Param.GetScalar <float> ("SimulationControl.MaximumTopGridTimeStep");

    MetaData.TimeLastRestartDump    = Param.GetScalar <float> ("InternalParameters.outputLabeling.TimeLastRestartDump");
    MetaData.dtRestartDump          = Param.GetScalar <float> ("OutputControlParameters.restartDump.dtRestartDump");
    MetaData.TimeLastDataDump       = Param.GetScalar <FLOAT> ("InternalParameters.outputLabeling.TimeLastDataDump");
    MetaData.dtDataDump             = Param.GetScalar <FLOAT> ("OutputControlParameters.dataDump.dtDataDump");
    MetaData.TimeLastHistoryDump    = Param.GetScalar <FLOAT> ("InternalParameters.outputLabeling.TimeLastHistoryDump");
    MetaData.dtHistoryDump          = Param.GetScalar <FLOAT> ("dtHistoryDump");
 
    TracerParticleOn                      = Param.GetScalar <int> ("SimulationControl.TracerParticleOn"); // should be bool
    ParticleTypeInFile                    = Param.GetScalar <int> ("OutputControlParameters.ParticleTypeInFile"); // should be bool
    OutputParticleTypeGrouping            = Param.GetScalar <int> ("OutputControlParameters.OutputParticleTypeGrouping");
    MetaData.TimeLastTracerParticleDump   = Param.GetScalar <FLOAT> ("TimeLastTracerParticleDump"); //doc?
    MetaData.dtTracerParticleDump         = Param.GetScalar <FLOAT> ("dtTracerParticleDump");//doc?
    MetaData.TimeLastInterpolatedDataDump = Param.GetScalar <FLOAT> ("TimeLastInterpolatedDataDump");//doc?
    MetaData.dtInterpolatedDataDump       = Param.GetScalar <FLOAT> ("dtInterpolatedDataDump");//doc?

 
    Param.GetArray <FLOAT> ("OutputControlParameters.movieDump.NewMovieLeftEdge", MetaData.NewMovieLeftEdge);
    Param.GetArray <FLOAT> ("OutputControlParameters.movieDump.NewMovieLeftEdge", MetaData.NewMovieRightEdge);

    MetaData.CycleLastRestartDump    = Param.GetScalar <int> ("CycleLastRestartDump"); //not used
    MetaData.CycleSkipRestartDump    = Param.GetScalar <int> ("CycleSkipRestartDump"); //not used
    MetaData.CycleLastDataDump       = Param.GetScalar <int> ("OutputControlParameters.outputLabeling.CycleLastDataDump");
    MetaData.CycleSkipDataDump       = Param.GetScalar <int> ("OutputControlParameters.cycleDump.CycleSkipDataDump");
    MetaData.CycleLastHistoryDump    = Param.GetScalar <int> ("OutputControlParameters.outputLabeling.CycleLastHistoryDump");
    MetaData.CycleSkipHistoryDump    = Param.GetScalar <int> ("OutputControlParameters.cycleDump.CycleSkipHistoryDump");
    MetaData.CycleSkipGlobalDataDump = Param.GetScalar <int> ("CycleSkipGlobalDataDump");
    MetaData.OutputFirstTimeAtLevel  = Param.GetScalar <int> ("OutputControlParameters.outputTriggers.OutputFirstTimeAtLevel");
    MetaData.StopFirstTimeAtLevel    = Param.GetScalar <int> ("SimulationControl.StopFirstTimeAtLevel");

 
    /* Maximum density directed output */
    OutputOnDensity        = Param.GetScalar <int> ("OutputControlParameters.outputTriggers.OutputOnDensity");  // should be bool
    StartDensityOutputs    = Param.GetScalar <float> ("OutputControlParameters.outputTriggers.StartDensityOutputs");
    CurrentDensityOutput   = Param.GetScalar <float> ("OutputControlParameters.outputTriggers.CurrentDensityOutput");
    IncrementDensityOutput = Param.GetScalar <float> ("OutputControlParameters.outputTriggers.IncrementDensityOutput");


    /* Subcycle directed output */
    MetaData.SubcycleSkipDataDump = Param.GetScalar <int> ("OutputControlParameters.cycleDump.SubcycleSkipDataDump");
    MetaData.SubcycleLastDataDump = Param.GetScalar <int> ("InternalParameters.outputLabeling.SubcycleLastDataDump");
    MetaData.SubcycleNumber = Param.GetScalar <int> ("InternalParameters.SubcycleNumber");


    FileDirectedOutput = Param.GetScalar <int> ("OutputControlParameters.FileDirectedOutput");
    WriteBinaryHierarchy = Param.GetScalar <int> ("OutputControlParameters.WriteBinaryHierarchy");


    MetaData.RestartDumpNumber = Param.GetScalar <int> ("InternalParameters.outputLabeling.RestartDumpNumber");
    MetaData.DataDumpNumber = Param.GetScalar <int> ("InternalParameters.outputLabeling.DataDumpNumber");
    MetaData.HistoryDumpNumber = Param.GetScalar <int> ("InternalParameters.outputLabeling.HistoryDumpNumber");
    MetaData.TracerParticleDumpNumber = Param.GetScalar <int> ("InternalParameters.outputLabeling.TracerParticleDumpNumber");
 
    if (sscanf(line, "RestartDumpName      = %s", dummy) == 1) MetaData.RestartDumpName = dummy;
    if (sscanf(line, "DataDumpName         = %s", dummy) == 1) MetaData.DataDumpName = dummy;
    if (sscanf(line, "HistoryDumpName      = %s", dummy) == 1) MetaData.HistoryDumpName = dummy;
    if (sscanf(line, "TracerParticleDumpName = %s", dummy) == 1) MetaData.TracerParticleDumpName = dummy;
    if (sscanf(line, "RedshiftDumpName     = %s", dummy) == 1) MetaData.RedshiftDumpName = dummy;
 
    if (sscanf(line, "RestartDumpDir      = %s", dummy) == 1) MetaData.RestartDumpDir = dummy;
    if (sscanf(line, "DataDumpDir         = %s", dummy) == 1) MetaData.DataDumpDir = dummy;
    if (sscanf(line, "HistoryDumpDir      = %s", dummy) == 1) MetaData.HistoryDumpDir = dummy;
    if (sscanf(line, "TracerParticleDumpDir = %s", dummy) == 1) MetaData.TracerParticleDumpDir = dummy;
    if (sscanf(line, "RedshiftDumpDir     = %s", dummy) == 1) MetaData.RedshiftDumpDir = dummy;
 
    if (sscanf(line, "LocalDir            = %s", dummy) == 1) MetaData.LocalDir = dummy;
    if (sscanf(line, "GlobalDir           = %s", dummy) == 1) MetaData.GlobalDir = dummy;

    LoadBalancing = Param.GetScalar <int> ("SimulationControl.optimization.LoadBalancing");
    ResetLoadBalancing = Param.GetScalar <int> ("SimulationControl.optimization.ResetLoadBalancing");
    LoadBalancingCycleSkip = Param.GetScalar <int> ("SimulationControl.optimization.LoadBalancingCycleSkip");
    LoadBalancingMinLevel = Param.GetScalar <int> ("SimulationControl.optimization.LoadBalancingMinLevel");
    LoadBalancingMaxLevel = Param.GetScalar <int> ("SimulationControl.optimization.LoadBalancingMaxLevel");
 

#ifdef UNUSED // need to re-visit
    if (sscanf(line, "TimeActionType[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2) {
      ret++; TimeActionType[dim] = int_dummy;
      if (dim >= MAX_TIME_ACTIONS-1) {
	ENZO_VFAIL("Time action %"ISYM" > maximum allowed.\n", dim)
      }
    }
    if (sscanf(line, "TimeActionRedshift[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionRedshift[%"ISYM"] = %"PSYM, &dim,
		    TimeActionRedshift+dim);
    if (sscanf(line, "TimeActionTime[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionTime[%"ISYM"] = %"PSYM, &dim,
		    TimeActionTime+dim);
    if (sscanf(line, "TimeActionParameter[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionParameter[%"ISYM"] = %"FSYM, &dim,
		    TimeActionParameter+dim);
#endif 

    MetaData.StaticHierarchy = Param.GetScalar <int> ("SimulationControl.amr.StaticHierarchy");  // should be bool
 
    MetaData.TopGridRank = Param.GetScalar <int> ("SimulationControl.domain.TopGridRank");
    Param.GetArray <int> ("TopGridDimensions", MetaData.TopGridDims);
 
    MetaData.GravityBoundary = Param.GetScalar <int> ("PhysicsParameters.TopGridGravityBoundary");
 
#ifdef TRANSFER
    if (sscanf(line, "RadHydroParamfile = %s", dummy) == 1) MetaData.RadHydroParameterFname = dummy;
#endif
    ImplicitProblem = Param.GetScalar <int> ("PhysicsParameters.radiationTransfer.ImplicitProblem"); // should be bool
    RadiativeTransferFLD = Param.GetScalar <int> ("PhysicsParameters.radiationTransfer.RadiativeTransferFLD");
#ifdef EMISSIVITY
    StarMakerEmissivityField = Param.GetScalar <int> ("StarMakerEmissivityField");
    uv_param = Param.GetScalar <float> ("uv_param");
#endif

    MetaData.ParticleBoundaryType = Param.GetScalar <int> ("PhysicsParameters.ParticleBoundaryType");
    MetaData.NumberOfParticles = Param.GetScalar <int> ("InternalParameters.NumberOfParticles");
 
    MetaData.CourantSafetyNumber = Param.GetScalar <float> ("PhysicsParameters.hydro.CourantSafetyNumber");
    MetaData.PPMFlatteningParameter = Param.GetScalar <int> ("PhysicsParameters.hydro.PPMFlatteningParameter");  // should be bool
    MetaData.PPMDiffusionParameter = Param.GetScalar <int> ("PhysicsParameters.hydro.PPMDiffusionParameter");  // should be bool
    MetaData.PPMSteepeningParameter = Param.GetScalar <int> ("PhysicsParameters.hydro.PPMSteepeningParameter");  // should be bool

// **********************************************************************
 
    /* read global Parameters */


    ProblemType = Param.GetScalar <int> ("ProblemType");

#ifdef NEW_PROBLEM_TYPES
    if (sscanf(line, "ProblemTypeName = %s", dummy) == 1) {
      ProblemTypeName = dummy;
      ProblemType = -978;
      ret = 1;
    }
#endif
   
    HydroMethod=Param.GetScalar <int> ("HydroMethod");

    if (HydroMethod==MHD_RK) useMHD = 1;

 
    huge_number=Param.GetScalar <float> ("huge_number");
    tiny_number=Param.GetScalar <float> ("tiny_number");
    Gamma=param.GetScalar <float> ("Gamma");
    PressureFree                           =Param.GetScalar <float> ("PressureFree");
    RefineBy                               =Param.GetScalar <float> ("RefineBy");
    MaximumRefinementLevel                 =Param.GetScalar <int> ("MaximumRefinementLevel");
    MaximumGravityRefinementLevel          =Param.GetScalar <int> ("MaximumGravityRefinementLevel");
    MaximumParticleRefinementLevel         =Param.GetScalar <int> ("MaximumParticleRefinementLevel");
    Param.GetArray <int> ("CellFlaggingMethod",CellFlaggingMethod+0,CellFlaggingMethod+1,
			  CellFlaggingMethod+2,CellFlaggingMethod+3,CellFlaggingMethod+4,
			  CellFlaggingMethod+5,CellFlaggingMethod+6);

    FluxCorection                          =Param.GetScalar <int> ("Fluxcorrection");
    InterpolationMethod                    =Param.GetScalar <int> ("InterpolationMethod");
    ConservativeInterpolation              =Param.GetScalar <int> ("ConservativeInterpolation");
    MinimumEfficiency                      =Param.GetScalar <float> ("MinimumEfficiency");
    SubgridSizeAutoAdjust                  =Param.GetScalar <int> ("SubgridSizeAutoAdjust");
    OptimalSubgridsPerProcessor            =Param.GetScalar <int> ("OptimalSubgridsPerProcessor");
    MinimumSubgridEdge                     =Param.GetScalar <int> ("MinimumSubgridEdge");
    MaximumSubgridSize                     =Param.GetScalar <int> ("MaximumSubgridSize");
    NumberOfBufferZones                    =Param.GetScalar <int> ("NumberOfBufferZones");
    FastSiblingLocatorEntireDomain         =Param.GetScalar <int> ("FastSiblingLocatorEntireDomain");
    MustRefineRegionMinRefinementLevel     =Param.GetScalar <int> ("MustRefineRegionMinRefinementLevel");
    MetallicityRefinementMinLevel          =Param.GetScalar <int> ("MetallicityRefinementMinLevel");  
    MetallicityRefinementMinMetallicity    =Param.GetScalar <int> ("MetallicityRefinementMinMetallicity"); 
    MetallicityRefinementMinDensity        =Param.GetScalar <float> ("MetallicityRefinementMinDensity"); 
    Param.GetArray <FLOAT> ("DomainLeftEdge",DomainLeftEdge,DomainLeftEdge+1,DomainLeftEdge+2);
    Param.GetArray <FLOAT> ("DomainRighttEdge",DomainRightEdge,DomainRightEdge+1,DomainRightEdge+2);
    Param.GetArray <float> ("GridVelocity",GridVelocity,GridVelocity+1,GridVelocity+2);

   RefineRegionAutoAdjust                  =Param.GetScalar <int> ("RefineRegionAutoAdjust"); 
   Param.GetArray <FLOAT> ("RefineRegionLeftEdge",RefineRegionLeftEdge,RefineRegionLeftEdge+1,RefineRegionLeftEdge+2);
   Param.GetArray <FLOAT> ("RefineRegionRightEdge",RefineRegionRightEdge,RefineRegionRightEdge+1,RefineRegionRightEdge+2);
   Param.GetArray <FLOAT> ("MustRefineRegionLeftEdge",MustRefineRegionLeftEdge,MustRefineRegionLeftEdge+1,MustRefineRegionLeftEdge+2);
   Param.GetArray <FLOAT> ("MustRefineRegionRightEdge",MustRefineRegionRightEdge,MustRefineRegionRightEdge+1,MustRefineRegionRightEdge+2);

 
    
// **********************************************************************

    /* Read evolving RefineRegion */

    RefineRegionTimeType = Param.GetScalar <int> ("SimulationControl.amr.RefineRegionTimeType");
    if (sscanf(line, "RefineRegionFile = %s", dummy) == 1) {
      RefineRegionFile = dummy;
      ret++;
    }
    
    if (sscanf(line, "DataLabel[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataLabel[dim] = dummy;
    if (sscanf(line, "DataUnits[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataUnits[dim] = dummy;
 
    UniformGravity          = Param.GetScalar <int> ("PhysicsParameters.gravity.UniformGravity");
    UniformGravityDirection = Param.GetScalar <int> ("PhysicsParameters.gravity.UniformGravityDirection");
    UniformGravityConstant  = Param.GetScalar <float> ("PhysicsParameters.gravity.UniformGravityConstant");
 
    PointSourceGravity           = Param.GetScalar <int> ("PhysicsParameters.gravity.PointSourceGravity"); // should be bool
    Param.GetArray <FLOAT> ("PointSourceGravityPosition", PointSourceGravityPosition); // why is the syntax like this?
    PointSourceGravityConstant   = Param.GetScalar <float> ("PhysicsParameters.gravity.PointSourceGravityConstant");
    PointSourceGravityCoreRadius = Param.GetScalar <float> ("PhysicsParameters.gravity.PointSourceGravityCoreRadius");
 
    ExternalGravity              = Param.GetScalar <int> ("PhysicsParameters.gravity.ExternalGravity");

    SelfGravity                  = Param.GetScalar <int> ("PhysicsParameters.gravity.SelfGravity");
    SelfGravityGasOff            = Param.GetScalar <int> ("PhysicsParameters.gravity.SelfGravityGasOff");
    AccretionKernel              = Param.GetScalar <int> ("PhysicsParameters.AccretionKernel");
    GravitationalConstant        = Param.GetScalar <float> ("PhysicsParameters.gravity.GravitationalConstant");

    S2ParticleSize                 = Param.GetScalar <float> ("S2ParticleSize");
    GravityResolution              = Param.GetScalar <float> ("PhysicsParameters.gravity.GravityResolution");
    ComputePotential               = Param.GetScalar <int> ("PhysicsParameters.gravity.ComputePotential"); // should be bool
    PotentialIterations            = Param.GetScalar <int> ("PhysicsParameters.gravity.PotentialIterations");
    WritePotential                 = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.WritePotential");  // should be bool
    BaryonSelfGravityApproximation = Param.GetScalar <int> ("PhysicsParameters.gravity.BaryonSelfGravityApproximation");  // should be bool
 
    GreensFunctionMaxNumber = Param.GetScalar <int> ("PhysicsParameters.gravity.GreensFunctionMaxNumber");
    GreensFunctionMaxSize   = Param.GetScalar <int> ("PhysicsParameters.gravity.GreensFunctionMaxSize");
 
    DualEnergyFormalism                              = Param.GetScalar <int> ("PhysicsParameters.hydro.DualEnergyFormalism"); // should be bool
    DualEnergyFormalismEta1                          = Param.GetScalar <float> ("PhysicsParameters.hydro.DualEnergyFormalismEta1");
    DualEnergyFormalismEta2                          = Param.GetScalar <float> ("PhysicsParameters.hydro.DualEnergyFormalismEta2");
    ParticleCourantSafetyNumber                      = Param.GetScalar <float> ("PhysicsParameters.hydro.ParticleCourantSafetyNumber");
    RootGridCourantSafetyNumber                      = Param.GetScalar <int> ("PhysicsParameters.hydro.RootGridCourantSafetyNumber");
    RandomForcing                                    = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.RandomForcing");
    RandomForcingEdot                                = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.RandomForcingEdot");
    RandomForcingMachNumber                          = Param.GetScalar <float> ("RandomForcingMachNumber");
    RadiativeCooling                                 = Param.GetScalar <int> ("RadiativeCooling"); // should be bool
    RadiativeCoolingModel                            = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.RadiativeCoolingModel");
    GadgetEquilibriumCooling                         = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.GadgetEquilibriumCooling");
    MultiSpecies                                     = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.MultiSpecies");
    PrimordialChemistrySolver                        = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.PrimordialChemistrySolver");
    CIECooling                                       = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.CIECooling"); // should be bool
    H2OpticalDepthApproximation                      = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.H2OpticalDepthApproximation"); // should be bool
    ThreeBodyRate = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.ThreeBodyRate"); // should be bool

    if (sscanf(line, "CloudyCoolingGridFile = %s", dummy) == 1) {
      CloudyCoolingData.CloudyCoolingGridFile = dummy;
      ret++;
    }
    CloudyCoolingData.IncludeCloudyHeating           = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.cloudyCooling.IncludeCloudyHeating"); // should be bool
    CloudyCoolingData.IncludeCloudyMMW               = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.cloudyCooling.IncludeCloudyMMW"); // should be bool
    CloudyCoolingData.CMBTemperatureFloor            = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.cloudyCooling.CMBTemperatureFloor"); // should be bool
    CloudyCoolingData.CloudyMetallicityNormalization = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.cloudyCooling.CloudyMetallicityNormalization");
    CloudyCoolingData.CloudyElectronFractionFactor   = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.cloudyCooling.CloudyElectronFractionFactor");
    MetalCooling = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.MetalCooling"); // should be bool
    if (sscanf(line, "MetalCoolingTable = %s", dummy) == 1) {
      MetalCoolingTable = dummy;
      ret++;
    }

    CRModel = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.CRModel");
    ShockMethod = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.ShockMethod");
    ShockTemperatureFloor = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.ShockTemperatureFloor");
    StorePreShockFields = Param.GetScalar <int> ("PhysicsParameters.miscellaneous.StorePreShockFields");

    RadiationFieldType = Param.GetScalar <int> ("PhysicsParameters.radiationField.RadiationFieldType");
    TabulatedLWBackground = Param.GetScalar <int> ("TabulatedLWBackground"); // should be bool
    AdjustUVBackground = Param.GetScalar <int> ("PhysicsParameters.radiationField.AdjustUVBackground");
    SetUVBAmplitude = Param.GetScalar <float> ("PhysicsParameters.radiationField.SetUVBAmplitude");
    SetHeIIHeatingScale = Param.GetScalar <float> ("PhysicsParameters.radiationField.SetHeIIHeatingScale");
    RadiationFieldLevelRecompute   = Param.GetScalar <int> ("PhysicsParameters.radiationField.RadiationFieldLevelRecompute"); // should be bool
    RadiationData.RadiationShield = Param.GetScalar <int> ("PhysicsParameters.radiationField.RadiationShield"); // should be bool
    CoolData.f3 = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationSpectrumNormalization");
    CoolData.alpha0 = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationSpectrumSlope");
    CoolData.f0to3 = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.CoolDataf0to3");
    CoolData.RadiationRedshiftOn = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationRedshiftOn");
    CoolData.RadiationRedshiftOff = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationRedshiftOff");
    CoolData.RadiationRedshiftFullOn = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationRedshiftFullOn");
    CoolData.RadiationRedshiftDropOff = Param.GetScalar <float> ("PhysicsParameters.radiationField.RadiationRedshiftDropOff");
    CoolData.HydrogenFractionByMass = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.HydrogenFractionByMass");
    CoolData.DeuteriumToHydrogenRatio = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.DeuteriumToHydrogenRatio");
    CoolData.NumberOfTemperatureBins = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.NumberOfTemperatureBins");
    CoolData.ih2co = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.CoolDataIh2co"); // should be bool
    CoolData.ipihdt = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.CoolDataIpiht"); // should be bool
    CoolData.TemperatureStart = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.TemperatureStart");
    CoolData.TemperatureEnd = Param.GetScalar <float> ("PhysicsParameters.atomicPhysics.TemperatureEnd");
    CoolData.comp_xray = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.CoolDataCompXray"); // should be bool
    CoolData.temp_xray = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.CoolDataTempXray"); // should be bool
    RateData.CaseBRecombination = Param.GetScalar <int> ("RateDataCaseBRecombination"); // should be bool
    PhotoelectricHeating = Param.GetScalar <int> ("PhysicsParameters.atomicPhysics.PhotoelectricHeating"); // should be bool
    
    OutputCoolingTime           = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.OutputCoolingTime"); // should be bool 
    OutputTemperature           = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.OutputTemperature"); // should be bool 
    OutputSmoothedDarkMatter    = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.OutputSmoothedDarkMatter"); // should be bool 
    SmoothedDarkMatterNeighbors = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.SmoothedDarkMatterNeighbors"); // should be bool 
    OutputGriddedStarParticle   = Param.GetScalar <int> ("OutputControlParameters.supplementalFields.OutputGriddedStarParticle"); // should be bool 

    ZEUSQuadraticArtificialViscosity = Param.GetScalar <float> ("PhysicsParameters.hydro.ZEUSQuadraticArtificialViscosity");
    ZEUSLinearArtificialViscosity     = Param.GetScalar <float> ("PhysicsParameters.hydro.ZEUSLinearArtificialViscosity");
 
    UseMinimumPressureSupport           = Param.GetScalar <int> ("PhysicsParameters.hydro.UseMinimumPressureSupport"); // should be bool 
    MinimumPressureSupportParameter     = Param.GetScalar <float> ("PhysicsParameters.hydro.MinimumPressureSupportParameter");
    RefineByJeansLengthSafetyFactor     = Param.GetScalar <float> ("SimulationControl.amr.RefineByJeansLengthSafetyFactor");
    JeansRefinementColdTemperature      = Param.GetScalar <float> ("SimulationControl.amr.JeansRefinementColdTemperature");
    RefineByResistiveLengthSafetyFactor = Param.GetScalar <float> ("SimulationControl.amr.RefineByResistiveLengthSafetyFactor");
    MustRefineParticlesRefineToLevel    = Param.GetScalar <int> ("SimulationControl.amr.MustRefineParticlesRefineToLevel");
    MustRefineParticlesRefineToLevelAutoAdjust = Param.GetScalar <int> ("SimulationControl.amr.MustRefineParticlesRefineToLevelAutoAdjust"); // should be bool 
    MustRefineParticlesMinimumMass      = Param.GetScalar <float> ("SimulationControl.amr.MustRefineParticlesMinimumMass");
    ParticleTypeInFile                  = Param.GetScalar <int> ("OutputControlParameters.ParticleTypeInFile"); // should be bool 

    int NumberOfStaticRefineRegions = Param.Size("StaticRefineRegion.Regions");
    if (NumberOfStaticRefineRegions > MAX_STATIC_REGIONS-1) {
      ENZO_VFAIL("You've exceeded the maximum number of StaticRefineRegions (%d)!\n",MAX_STATIC_REGIONS)
    }

    std::string StaticRefineRegionNames[NumberOfStaticRefineRegions];
    Param.GetArray <std::string> ("StaticRefineRegion.Regions",StaticRefineRegionNames);

    for (i = 0; i < NumberOfStaticRefineRegions; i++) {
      StaticRefineRegionLevel[i] = Param.GetScalar <int> ("StaticRefineRegion.%s.Level",StaticRefineRegionNames[i]);
      Param.GetArray <FLOAT> ("StaticRefineRegion.%s.LeftEdge",StaticRefineRegionLeftEdge[i]);
      Param.GetArray <FLOAT> ("StaticRefineRegion.%s.RightEdge",StaticRefineRegionRightEdge[i]);
    }
    
    ParallelRootGridIO = Param.GetScalar <int> ("ParallelRootGridIO"); // should be bool
    ParallelParticleIO = Param.GetScalar <int> ("ParallelParticleIO"); // should be bool
 
    Unigrid          = Param.GetScalar <int> ("Unigrid"); // should be bool
    UnigridTranspose = Param.GetScalar <int> ("UnigridTranspose"); // should be bool
    NumberOfRootGridTilesPerDimensionPerProcessor = Param.GetScalar <int> ("NumberOfRootGridTilesPerDimensionPerProcessor");
 
    PartitionNestedGrids = Param.GetScalar <int> ("PartitionNestedGrids");
    StaticPartitionNestedGrids = Param.GetScalar <int> ("StaticPartitionNestedGrids"); // should be bool
 
    ExtractFieldsOnly = Param.GetScalar <int> ("ExtractFieldsOnly"); // should be bool
 
    debug1 = Param.GetScalar <int> ("Debug1"); // should be bool

    debug2 = Param.GetScalar <int> ("Debug2"); // should be bool

    MemoryLimit = Param.GetScalar <int> ("MemoryLimit");

#ifdef STAGE_INPUT
    StageInput = Param.GetScalar <int> ("StageInput");
#endif

#ifdef OOC_BOUNDARY

    ExternalBoundaryIO = Param.GetScalar <int> ("ExternalBoundaryIO"); // should be bool

    ExternalBoundaryTypeIO = Param.GetScalar <int> ("ExternalBoundaryTypeIO"); // should be bool

    ExternalBoundaryValueIO = Param.GetScalar <int> ("ExternalBoundaryValueIO");
 // should be bool

    SimpleConstantBoundary = Param.GetScalar <int> ("SimpleConstantBoundary"); // should be bool

#endif

    Param.GetArray <int> ("SlopeFlaggingFields", SlopeFlaggingFields);

    Param.GetArray <float> ("MinimumSlopeForRefinement", MinimumSlopeForRefinement);
 
    Param.GetArray <float> ("MinimumOverDensityForRefinement", MinimumOverDensityForRefinement);
 
    Param.GetArray <float> ("MinimumMassForRefinement", MinimumMassForRefinement);

    Param.GetArray <float> ("MinimumMassForRefinementLevelExponent", MinimumMassForRefinementLevelExponent);


    MinimumPressureJumpForRefinement = Param.GetScalar <float> ("MinimumPressureJumpForRefinement");
    MinimumShearForRefinement        = Param.GetScalar <float> ("MinimumShearForRefinement");
    MinimumEnergyRatioForRefinement  = Param.GetScalar <float> ("MinimumEnergyRatioForRefinement");
    ShockwaveRefinementMinMach       = Param.GetScalar <float> ("ShockwaveRefinementMinMach");
    ShockwaveRefinementMinVelocity   = Param.GetScalar <float> ("ShockwaveRefinementMinVelocity");
    ShockwaveRefinementMaxLevel      = Param.GetScalar <float> ("ShockwaveRefinementMaxLevel");
    ComovingCoordinates              = Param.GetScalar <int> ("ComovingCoordinates"); // should be bool
    StarParticleCreation             = Param.GetScalar <int> ("StarParticleCreation");
    BigStarFormation                 = Param.GetScalar <int> ("BigStarFormation"); // should be bool
    BigStarFormationDone             = Param.GetScalar <int> ("BigStarFormationDone"); // should be bool
    BigStarSeparation                = Param.GetScalar <float> ("BigStarSeparation");
    SimpleQ                          = Param.GetScalar <float> ("SimpleQ");
    SimpleRampTime                   = Param.GetScalar <float> ("SimpleRampTime");
    StarParticleFeedback             = Param.GetScalar <int> ("StarParticleFeedback");
    NumberOfParticleAttributes       = Param.GetScalar <int> ("NumberOfParticleAttributes");
    AddParticleAttributes            = Param.GetScalar <int> ("AddParticleAttributes"); // should be bool


    /* read data which defines the boundary conditions */
    Param.GetArray <int> ("LeftFaceBoundaryCondition", MetaData.LeftFaceBoundaryCondition);
    Param.GetArray <int> ("RightFaceBoundaryCondition", MetaData.RightFaceBoundaryCondition);
  
    if (sscanf(line, "BoundaryConditionName         = %s", dummy) == 1)
      MetaData.BoundaryConditionName = dummy;

    if (sscanf(line, "MetaDataIdentifier = %s", dummy) == 1) {
      MetaData.MetaDataIdentifier = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataSimulationUUID = %s", dummy) == 1) {
      MetaData.SimulationUUID = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataDatasetUUID = %s", dummy) == 1) {
      MetaData.RestartDatasetUUID = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataInitialConditionsUUID = %s", dummy) == 1) {
      MetaData.InitialConditionsUUID = dummy;
      ret++;
    }

    /* Check version number. */
 
    if (sscanf(line, "VersionNumber = %"FSYM, &TempFloat) == 1) {
      ret++;
      if (fabs(TempFloat - VERSION) >= 1.0e-3 &&
	  MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Warning: Incorrect version number.\n");
    }
 
    /* Read star particle parameters. */

    StarMakerOverDensityThreshold = Param.GetScalar <float> ("StarMakerOverDensityThreshold");
    StarMakerSHDensityThreshold      = Param.GetScalar <float> ("StarMakerSHDensityThreshold");
    StarMakerMassEfficiency     = Param.GetScalar <float> ("StarMakerMassEfficiency");
    StarMakerMinimumMass     = Param.GetScalar <float> ("StarMakerMinimumMass");
    StarMakerMinimumDynamicalTime     = Param.GetScalar <float> ("StarMakerMinimumDynamicalTime");
    StarMassEjectionFraction          = Param.GetScalar <float> ("StarMassEjectionFraction");
    StarMetalYield                    = Param.GetScalar <float> ("StarMetalYield");
    StarEnergyToThermalFeedback       = Param.GetScalar <float> ("StarEnergyToThermalFeedback");
    StarEnergyToStellarUV = Param.GetScalar <float> ("StarEnergyToStellarUV");
    StarEnergyToQuasarUV= Param.GetScalar <float> ("StarEnergyToQuasarUV");

= Param.GetScalar <int> ("");
= Param.GetScalar <int> ("");
= Param.GetScalar <int> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <int> ("");
= Param.GetScalar <int> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <float> ("");
= Param.GetScalar <float> ("");

= Param.GetScalar <> ("");

    ret += sscanf(line, "StarFeedbackDistRadius = %"ISYM, &StarFeedbackDistRadius);
    ret += sscanf(line, "StarFeedbackDistCellStep = %"ISYM, &StarFeedbackDistCellStep);

    ret += sscanf(line, "StarClusterUseMetalField = %"ISYM, 		  &StarClusterUseMetalField);
    ret += sscanf(line, "StarClusterMinDynamicalTime = %"FSYM, 		  &StarClusterMinDynamicalTime);
    ret += sscanf(line, "StarClusterHeliumIonization = %"ISYM, 		  &StarClusterHeliumIonization);
    ret += sscanf(line, "StarClusterUnresolvedModel = %"ISYM, 		  &StarClusterUnresolvedModel);
    ret += sscanf(line, "StarClusterIonizingLuminosity = %lf", 		  &StarClusterIonizingLuminosity);
    ret += sscanf(line, "StarClusterSNEnergy = %lf", &StarClusterSNEnergy);
    ret += sscanf(line, "StarClusterSNRadius = %"FSYM, &StarClusterSNRadius);
    ret += sscanf(line, "StarClusterFormEfficiency = %"FSYM, 		  &StarClusterFormEfficiency);
    ret += sscanf(line, "StarClusterMinimumMass = %"FSYM, 		  &StarClusterMinimumMass);
    ret += sscanf(line, "StarClusterCombineRadius = %"FSYM,		  &StarClusterCombineRadius);
    ret += sscanf(line, "StarClusterRegionLeftEdge = %"FSYM" %"FSYM" %"FSYM,
		  StarClusterRegionLeftEdge, StarClusterRegionLeftEdge+1, 
		  StarClusterRegionLeftEdge+2);
    ret += sscanf(line, "StarClusterRegionRightEdge = %"FSYM" %"FSYM" %"FSYM,
		  StarClusterRegionRightEdge, StarClusterRegionRightEdge+1, 
		  StarClusterRegionRightEdge+2);

    ret += sscanf(line, "PopIIIStarMass = %"FSYM, &PopIIIStarMass);
    ret += sscanf(line, "PopIIIInitialMassFunction = %"ISYM, 		  &PopIIIInitialMassFunction);
    ret += sscanf(line, "PopIIIInitialMassFunctionSeed = %"ISYM, 		  &PopIIIInitialMassFunctionSeed);
    ret += sscanf(line, "PopIIIInitialMassFunctionCalls = %"ISYM, 		  &PopIIIInitialMassFunctionCalls);
    ret += sscanf(line, "PopIIIMassRange = %"FSYM" %"FSYM,		  &PopIIILowerMassCutoff, &PopIIIUpperMassCutoff);
    ret += sscanf(line, "PopIIIInitialMassFunctionSlope = %"FSYM, 		  &PopIIIInitialMassFunctionSlope);
    ret += sscanf(line, "PopIIIBlackHoles = %"ISYM, &PopIIIBlackHoles);
    ret += sscanf(line, "PopIIIBHLuminosityEfficiency = %"FSYM, 		  &PopIIIBHLuminosityEfficiency);
    ret += sscanf(line, "PopIIIOverDensityThreshold = %"FSYM,		  &PopIIIOverDensityThreshold);
    ret += sscanf(line, "PopIIIH2CriticalFraction = %"FSYM,		  &PopIIIH2CriticalFraction);
    ret += sscanf(line, "PopIIIMetalCriticalFraction = %"FSYM,		  &PopIIIMetalCriticalFraction);
    ret += sscanf(line, "PopIIISupernovaRadius = %"FSYM, &PopIIISupernovaRadius);
    ret += sscanf(line, "PopIIISupernovaUseColour = %"ISYM, 		  &PopIIISupernovaUseColour);
    ret += sscanf(line, "PopIIISupernovaMustRefine = %"ISYM,		  &PopIIISupernovaMustRefine);
    ret += sscanf(line, "PopIIISupernovaMustRefineResolution = %"ISYM,		  &PopIIISupernovaMustRefineResolution);
    ret += sscanf(line, "PopIIIHeliumIonization = %"ISYM, 		  &PopIIIHeliumIonization);

    ret += sscanf(line, "PopIIIColorDensityThreshold = %"FSYM,		  &PopIIIColorDensityThreshold);
    ret += sscanf(line, "PopIIIColorMass = %"FSYM,		  &PopIIIColorMass);

    ret += sscanf(line, "MBHAccretion = %"ISYM, &MBHAccretion);
    ret += sscanf(line, "MBHAccretionRadius = %"FSYM, &MBHAccretionRadius);
    ret += sscanf(line, "MBHAccretingMassRatio = %"FSYM, &MBHAccretingMassRatio);
    ret += sscanf(line, "MBHAccretionFixedTemperature = %"FSYM, &MBHAccretionFixedTemperature);
    ret += sscanf(line, "MBHAccretionFixedRate = %"FSYM, &MBHAccretionFixedRate);
    ret += sscanf(line, "MBHTurnOffStarFormation = %"ISYM, &MBHTurnOffStarFormation);
    ret += sscanf(line, "MBHCombineRadius = %"FSYM, &MBHCombineRadius);
    ret += sscanf(line, "MBHMinDynamicalTime = %"FSYM, &MBHMinDynamicalTime);
    ret += sscanf(line, "MBHMinimumMass = %"FSYM, &MBHMinimumMass);

    ret += sscanf(line, "MBHFeedback = %"ISYM, &MBHFeedback);
    ret += sscanf(line, "MBHFeedbackRadiativeEfficiency = %"FSYM, &MBHFeedbackRadiativeEfficiency);
    ret += sscanf(line, "MBHFeedbackEnergyCoupling = %"FSYM, &MBHFeedbackEnergyCoupling);
    ret += sscanf(line, "MBHFeedbackMassEjectionFraction = %"FSYM, &MBHFeedbackMassEjectionFraction);
    ret += sscanf(line, "MBHFeedbackMetalYield = %"FSYM, &MBHFeedbackMetalYield);
    ret += sscanf(line, "MBHFeedbackThermalRadius = %"FSYM, &MBHFeedbackThermalRadius);
    ret += sscanf(line, "MBHFeedbackJetsThresholdMass = %"FSYM, &MBHFeedbackJetsThresholdMass);

    ret += sscanf(line, "MBHParticleIO = %"ISYM,
		  &MBHParticleIO);
    if (sscanf(line, "MBHParticleIOFilename = %s", dummy) == 1)
      MBHParticleIOFilename = dummy;
    if (sscanf(line, "MBHInsertLocationFilename = %s", dummy) == 1)
      MBHInsertLocationFilename = dummy;


    // **********************************************************************

    /* Read Movie Dump parameters */

    MovieSkipTimestep   = Param.GetScalar <int> ("MovieSkipTimestep");
    Movie3DVolumes      = Param.GetScalar <int> ("Movie3DVolumes");
    MovieVertexCentered = Param.GetScalar <int> ("MovieVertexCentered");
    NewMovieParticleOn  = Param.GetScalar <int> ("NewMovieParticleOn");
    Param.GetArray <int> ("MovieDataField", MovieDataField);
    NewMovieDumpNumber  = Param.GetScalar <int> ("NewMovieDumpNumber");
    if (sscanf(line, "NewMovieName = %s", dummy) == 1)
      NewMovieName = dummy;
    MetaData.MovieTimestepCounter = Param.GetScalar <int> ("MovieTimestepCounter");

    MultiMetals           = Param.GetScalar <int> ("MultiMetals");
    IsotropicConduction   = Param.GetScalar <int> ("IsotropicConduction"); // should be bool
    AnisotropicConduction = Param.GetScalar <int> ("AnisotropicConduction"); // should be bool
    IsotropicConductionSpitzerFraction   = Param.GetScalar <float> ("IsotropicConductionSpitzerFraction");
    AnisotropicConductionSpitzerFraction = Param.GetScalar <float> ("AnisotropicConductionSpitzerFraction");
    ConductionCourantSafetyNumber        = Param.GetScalar <float> ("ConductionCourantSafetyNumber");

    RadiativeTransfer           = Param.GetScalar <int> ("RadiativeTransfer"); // should be bool
    RadiationXRaySecondaryIon   = Param.GetScalar <int> ("RadiationXRaySecondaryIon"); // should be bool
    RadiationXRayComptonHeating = Param.GetScalar <int> ("RadiationXRayComptonHeating"); // should be bool


    /* Shearing Box Boundary parameters */

    AngularVelocity           = Param.GetScalar <float> ("AngularVelocity");
    VelocityGradient          = Param.GetScalar <float> ("VelocityGradient");
    ShearingVelocityDirection = Param.GetScalar <int> ("ShearingVelocityDirection");
    ShearingBoxProblemType    = Param.GetScalar <int> ("ShearingBoxProblemType");
    

#ifdef STAGE_INPUT
    sscanf(line, "LocalPath = %s\n", LocalPath);
    sscanf(line, "GlobalPath = %s\n", GlobalPath);
#endif

    /* Embedded Python */
    PythonTopGridSkip  = Param.GetScalar <int> ("PythonTopGridSkip"); // should be bool
    PythonSubcycleSkip = Param.GetScalar <int> ("PythonSubcycleSkip"); // should be bool

#ifdef USE_PYTHON
    NumberOfPythonCalls         = Param.GetScalar <int> ("NumberOfPythonCalls");
    NumberOfPythonTopGridCalls  = Param.GetScalar <int> ("NumberOfPythonTopGridCalls");
    NumberOfPythonSubcycleCalls = Param.GetScalar <int> ("NumberOfPythonSubcycleCalls");
#endif

    /* Inline halo finder */

    InlineHaloFinder             = Param.GetScalar <int> ("InlineHaloFinder"); // should be bool
    HaloFinderSubfind            = Param.GetScalar <int> ("HaloFinderSubfind"); // should be bool
    HaloFinderOutputParticleList = Param.GetScalar <int> ("HaloFinderOutputParticleList"); // should be bool
    HaloFinderRunAfterOutput     = Param.GetScalar <int> ("HaloFinderRunAfterOutput"); // should be bool
    HaloFinderLinkingLength      = Param.GetScalar <float> ("HaloFinderLinkingLength");
    HaloFinderMinimumSize        = Param.GetScalar <int> ("HaloFinderMinimumSize");
    HaloFinderCycleSkip          = Param.GetScalar <int> ("HaloFinderCycleSkip");
    HaloFinderTimestep           = Param.GetScalar <float> ("HaloFinderTimestep");
    HaloFinderLastTime           = Param.GetScalar <FLOAT> ("HaloFinderLastTime");

    /* This Block for Stanford Hydro */

    UseHydro = Param.GetScalar <int> ("UseHydro"); // should be bool


    /* Sink particles (for present day star formation) & winds */
    SinkMergeDistance       = Param.GetScalar <float> ("SinkMergeDistance");
    SinkMergeMass           = Param.GetScalar <float> ("SinkMergeMass");
    StellarWindFeedback     = Param.GetScalar <int> ("StellarWindFeedback"); // should be bool
    StellarWindTurnOnMass   = Param.GetScalar <float> ("StellarWindTurnOnMass");
    MSStellarWindTurnOnMass = Param.GetScalar <float> ("MSStellarWindTurnOnMass");

    VelAnyl = Param.GetScalar <int> ("VelAnyl"); // should be bool
    BAnyl   = Param.GetScalar <int> ("BAnyl"); // should be bool



    /* Read MHD Paramters */
    UseDivergenceCleaning            = Param.GetScalar <int> ("UseDivergenceCleaning"); // should be bool
    DivergenceCleaningBoundaryBuffer = Param.GetScalar <int> ("DivergenceCleaningBoundaryBuffer"); // should be bool
    DivergenceCleaningThreshold      = Param.GetScalar <float> ("DivergenceCleaningThreshold");
    PoissonApproximationThreshold    = Param.GetScalar <float> ("PoissonApproximationThreshold");
    PoissonBoundaryType              = Param.GetScalar <int> ("PoissonBoundaryType");
   

    AngularVelocity   = Param.GetScalar <float> ("AngularVelocity");
    VelocityGradient  = Param.GetScalar <float> ("VelocityGradient");
    UseDrivingField   = Param.GetScalar <int> ("UseDrivingField"); // should be bool
    DrivingEfficiency = Param.GetScalar <float> ("DrivingEfficiency");

    StringKick            = Param.GetScalar <float> ("StringKick");
    StringKickDimension   = Param.GetScalar <int> ("StringKickDimension");
    UsePhysicalUnit       = Param.GetScalar <int> ("UsePhysicalUnit"); // should be bool
    Theta_Limiter         = Param.GetScalar <float> ("Theta_Limiter");
    RKOrder               = Param.GetScalar <int> ("RKOrder");
    UseFloor              = Param.GetScalar <int> ("UseFloor"); // should be bool
    UseViscosity          = Param.GetScalar <int> ("UseViscosity"); // should be bool
    ViscosityCoefficient  = Param.GetScalar <float> ("ViscosityCoefficient");

    UseAmbipolarDiffusion = Param.GetScalar <int> ("UseAmbipolarDiffusion"); // should be bool
    UseResistivity        = Param.GetScalar <int> ("UseResistivity"); // should be bool
    SmallRho              = Param.GetScalar <float> ("SmallRho");
    SmallP                = Param.GetScalar <float> ("SmallP");
    SmallT                = Param.GetScalar <float> ("SmallT");
    MaximumAlvenSpeed     = Param.GetScalar <float> ("MaximumAlvenSpeed");
    Coordinate            = Param.GetScalar <int> ("Coordinate");
    RiemannSolver         = Param.GetScalar <int> ("RiemannSolver");
    RiemannSolverFallback = Param.GetScalar <int> ("RiemannSolverFallback"); // should be bool
    ConservativeReconstruction = Param.GetScalar <int> ("ConservativeReconstruction"); // should be bool
    PositiveReconstruction = Param.GetScalar <int> ("PositiveReconstruction"); // should be bool
    ReconstructionMethod  = Param.GetScalar <int> ("ReconstructionMethod");

    EOSType               = Param.GetScalar <int> ("EOSType");
    EOSSoundSpeed         = Param.GetScalar <float> ("EOSSoundSpeed");
    EOSCriticalDensity    = Param.GetScalar <float> ("EOSCriticalDensity");
    EOSGamma              = Param.GetScalar <float> ("EOSGamma");
    UseConstantAcceleration = Param.GetScalar <int> ("UseConstantAcceleration");
 // should be bool

    Param.GetArray <float> ("ConstantAcceleration", ConstantAcceleration);
    Mu = Param.GetScalar <float> ("Mu");
    DivBDampingLength = Param.GetScalar <float> ("DivBDampingLength");
    CoolingCutOffDensity1 = Param.GetScalar <float> ("CoolingCutOffDensity1");
    CoolingCutOffDensity2 = Param.GetScalar <float> ("CoolingCutOffDensity2");
    CoolingCutOffTemperature = Param.GetScalar <float> ("CoolingCutOffTemperature");
    CoolingPowerCutOffDensity1 = Param.GetScalar <float> ("CoolingPowerCutOffDensity1");
    CoolingPowerCutOffDensity2 = Param.GetScalar <float> ("CoolingPowerCutOffDensity2");
    UseH2OnDust = Param.GetScalar <int> ("UseH2OnDust"); // should be bool
    seCUDA = Param.GetScalar <int> ("UseCUDA"); // should be bool

    MoveParticlesBetweenSiblings = Param.GetScalar <int> ("MoveParticlesBetweenSiblings"); // should be bool
    ParticleSplitterIterations = Param.GetScalar <int> ("ParticleSplitterIterations");
    ParticleSplitterChildrenParticleSeparation = Param.GetScalar <float> ("ParticleSplitterChildrenParticleSeparation");
    ResetMagneticField = Param.GetScalar <int> ("ResetMagneticField"); // should be bool
    Param.GetArray <float> ("ResetMagneticFieldAmplitude", ResetMagneticFieldAmplitude);

    /* If the dummy char space was used, then make another. */
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
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

    if (strstr(line, "\"\"\"")              ) comment_count++;

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);
 
  }

  /* clean up */
 
  delete [] dummy;
  rewind(fptr);

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
    MustRefineParticlesMinimumMass /= POW(1/(float(MetaData.TopGridDims[0])
				       *POW(float(RefineBy), float(MustRefineParticlesRefineToLevel))),3);


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

    ret += sscanf(line, "IsothermalSoundSpeed = %"GSYM, &IsothermalSoundSpeed);
    ret += sscanf(line, "RefineByJeansLengthUnits = %"ISYM, &RefineByJeansLengthUnits);
 
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
