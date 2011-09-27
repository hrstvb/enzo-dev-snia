/***********************************************************************
/
/  SETS THE DEFAULT GLOBAL VALUES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       February 29th, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// Set default values for global variables that are not set in
// defaults.cfg.
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "StarParticleData.h"
 
/* character strings */
 
char DefaultDimUnits[] = "cm";
char *DefaultDimLabel[] = {"x", "y", "z"};
  
int SetDefaultGlobalValues(TopGridData &MetaData)
{
 
  /* declarations */
 
  const float Pi = 3.14159;
  int dim, i, j;
 
  /* set the default MetaData values. */
 
  MetaData.WroteData           = FALSE;
 
  MetaData.CycleSkipRestartDump = 0;
  MetaData.CycleSkipDataDump    = 0;
  MetaData.CycleSkipHistoryDump = 0;
  MetaData.CycleSkipGlobalDataDump = 0; //AK
  
  // MetaData.NumberOfOutputsBeforeExit = 0;
  // MetaData.OutputsLeftBeforeExit     = 0;


  CoresPerNode = 1;
  PreviousMaxTask = 0;

  // for (i = 0;i < MAX_DEPTH_OF_HIERARCHY;i++) {
  //   RebuildHierarchyCycleSkip[i] = 1;
  // }

 
  for (i = 0; i < MAX_CUBE_DUMPS; i++) {
    CubeDumps[i] = NULL;
  }
  
  MetaData.FirstTimestepAfterRestart = TRUE;
 
 
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DimLabels[dim]                  = DefaultDimLabel[dim];
    DimUnits[dim]                   = DefaultDimUnits;
    // ExternalGravityPosition[dim]    = 0.0;
    // ExternalGravityOrientation[dim] = 0.0;
  }
 
  // for (i = 0; i < MAX_STATIC_REGIONS; i++) {
  //   AvoidRefineRegionLevel[i]  = INT_UNDEFINED;
  //   for (dim = 0; dim < MAX_DIMENSION; dim++) {
  //     AvoidRefineRegionLeftEdge[i][dim] = FLOAT_UNDEFINED;
  //     AvoidRefineRegionRightEdge[i][dim] = FLOAT_UNDEFINED;
  //   }
  // }

  /* For evolving refinement regions. */
  for (i = 0; i < MAX_REFINE_REGIONS; i++) {
    EvolveRefineRegionTime[i] = FLOAT_UNDEFINED;
    for (j = 0; j < MAX_DIMENSION; j++) {
      EvolveRefineRegionLeftEdge[i][j]  = FLOAT_UNDEFINED;
      EvolveRefineRegionRightEdge[i][j] = FLOAT_UNDEFINED;
    }
  }

  //  DatabaseLocation = NULL;

 

  ExternalBoundaryIO          = FALSE;
  ExternalBoundaryTypeIO      = FALSE;
  ExternalBoundaryValueIO     = FALSE;
  SimpleConstantBoundary      = FALSE;

  debug1                      = 0;
  debug2                      = 0;

  CurrentMaximumDensity            = -999;
 
#ifdef STAGE_INPUT
  StageInput                  = 0;
#endif

  First_Pass                  = 0;

  // ExternalGravityConstant = 0.0;
  ExternalGravityDensity      = 0.0;
  ExternalGravityRadius       = 0.0;

  CopyGravPotential           = FALSE;             // off
 
  StarMakerEmissivityField    = 0;                 // off
  uv_param                    = 1.0e-5;            // mid-range value from Razoumov Norman 2002

  // NoMultiSpeciesButColors     = FALSE;             // off
  // H2FormationOnDust           = FALSE;
  GloverChemistryModel        = 0;                 // 0ff
  GloverRadiationBackground   = 0;
  GloverOpticalDepth          = 0;

  // RadiationFieldRedshift      = 0.0;
  // PhotoelectricHeatingRate    = 8.5e-26;           // ergs cm-3 s-1

  CoolData.alpha0             = 1.5;               // radiation spectral slope
  CoolData.f3                 = 1.0e-21;           // radiation normalization
  // CoolData.SolarMetalFractionByMass = 0.02041;
  // RateData.NumberOfDustTemperatureBins = 250;
  // RateData.DustTemperatureStart    = 1.0;
  // RateData.DustTemperatureEnd      = 1500.0;

  // OutputDustTemperature = FALSE;

 
  StarFeedbackDistCellStep         = 0;
  StarFeedbackDistTotalCells       = 1;

  IMFData                          = NULL;

  OutputWhenJetsHaveNotEjected     = FALSE;

  LastSupernovaTime                = FLOAT_UNDEFINED;

  ran1_init = 0;
  rand_init = 0;

  NSpecies		     = INT_UNDEFINED;
  NColor		     = INT_UNDEFINED;
  NEQ_HYDRO		     = 5;
  NEQ_MHD		     = 9;
  SmallEint		     = 1e-30;

  iden	= 0;
  ivx	= 1;
  ivy	= 2;
  ivz	= 3;
  ietot = 4;
  ieint = 0;
  iBx	= 5;
  iBy	= 6;
  iBz	= 7;
  iPhi	= 8;
  iD	= 0;
  iS1	= 1;
  iS2	= 2;
  iS3	= 3;
  iEtot = 4;
  iEint = 0;

  /* End of Stanford Hydro additions */

  /* test problem values */
  TestProblemData.HydrogenFractionByMass = 0.76;

  /* The DToHRatio is by mass in the code, so multiply by 2. */
  TestProblemData.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998 

  // multispecies default values assume completely neutral gas with primordial D/H ratio
  TestProblemData.MultiSpecies = 0;
  TestProblemData.HI_Fraction = 1.0;
  TestProblemData.HII_Fraction = tiny_number;
  TestProblemData.HeI_Fraction = 1.0;
  TestProblemData.HeII_Fraction = tiny_number;
  TestProblemData.HeIII_Fraction = tiny_number;
  TestProblemData.HM_Fraction = tiny_number;
  TestProblemData.H2I_Fraction = tiny_number;
  TestProblemData.H2II_Fraction = tiny_number;
  TestProblemData.DI_Fraction = 2.0*3.4e-5;
  TestProblemData.DII_Fraction = tiny_number;
  TestProblemData.HDI_Fraction = tiny_number;

  // This is for ionized gas (within a supernova blast, for example)
  TestProblemData.HI_Fraction_Inner = 1.0;
  TestProblemData.HII_Fraction_Inner = tiny_number;
  TestProblemData.HeI_Fraction_Inner = 1.0;
  TestProblemData.HeII_Fraction_Inner = tiny_number;
  TestProblemData.HeIII_Fraction_Inner = tiny_number;
  TestProblemData.HM_Fraction_Inner = tiny_number;
  TestProblemData.H2I_Fraction_Inner = tiny_number;
  TestProblemData.H2II_Fraction_Inner = tiny_number;
  TestProblemData.DI_Fraction_Inner = 2.0*3.4e-5;
  TestProblemData.DII_Fraction_Inner = tiny_number;
  TestProblemData.HDI_Fraction_Inner = tiny_number;

  TestProblemData.UseMetallicityField = 0;
  TestProblemData.MetallicityField_Fraction = tiny_number;

  TestProblemData.UseMassInjection = 0;
  TestProblemData.InitialHydrogenMass = tiny_number;
  TestProblemData.InitialDeuteriumMass = tiny_number;
  TestProblemData.InitialHeliumMass = tiny_number;
  TestProblemData.InitialMetalMass = tiny_number;

  TestProblemData.MultiMetals = 0;
  TestProblemData.MultiMetalsField1_Fraction = tiny_number;
  TestProblemData.MultiMetalsField2_Fraction = tiny_number;

  TestProblemData.GloverChemistryModel = 0;
  // This is for the gas in the surrounding medium, for the blast wave problem.
  TestProblemData.CI_Fraction = tiny_number;
  TestProblemData.CII_Fraction = tiny_number;
  TestProblemData.OI_Fraction = tiny_number;
  TestProblemData.OII_Fraction = tiny_number;
  TestProblemData.SiI_Fraction = tiny_number;
  TestProblemData.SiII_Fraction = tiny_number;
  TestProblemData.SiIII_Fraction = tiny_number;
  TestProblemData.CHI_Fraction = tiny_number;
  TestProblemData.CH2I_Fraction = tiny_number;
  TestProblemData.CH3II_Fraction = tiny_number;
  TestProblemData.C2I_Fraction = tiny_number;
  TestProblemData.COI_Fraction = tiny_number;
  TestProblemData.HCOII_Fraction = tiny_number;
  TestProblemData.OHI_Fraction = tiny_number;
  TestProblemData.H2OI_Fraction = tiny_number;
  TestProblemData.O2I_Fraction = tiny_number;

  // This is for the gas in the region where the blast wave energy is injected
  TestProblemData.CI_Fraction_Inner = tiny_number;
  TestProblemData.CII_Fraction_Inner = tiny_number;
  TestProblemData.OI_Fraction_Inner = tiny_number;
  TestProblemData.OII_Fraction_Inner = tiny_number;
  TestProblemData.SiI_Fraction_Inner = tiny_number;
  TestProblemData.SiII_Fraction_Inner = tiny_number;
  TestProblemData.SiIII_Fraction_Inner = tiny_number;
  TestProblemData.CHI_Fraction_Inner = tiny_number;
  TestProblemData.CH2I_Fraction_Inner = tiny_number;
  TestProblemData.CH3II_Fraction_Inner = tiny_number;
  TestProblemData.C2I_Fraction_Inner = tiny_number;
  TestProblemData.COI_Fraction_Inner = tiny_number;
  TestProblemData.HCOII_Fraction_Inner = tiny_number;
  TestProblemData.OHI_Fraction_Inner = tiny_number;
  TestProblemData.H2OI_Fraction_Inner = tiny_number;
  TestProblemData.O2I_Fraction_Inner = tiny_number;

  TestProblemData.MinimumHNumberDensity = 1;
  TestProblemData.MaximumHNumberDensity = 1e6;
  TestProblemData.MinimumMetallicity    = 1e-6;
  TestProblemData.MaximumMetallicity    = 1;
  TestProblemData.MinimumTemperature    = 10;
  TestProblemData.MaximumTemperature    = 1e7;
  TestProblemData.ResetEnergies         = 1;

  // This should only be false for analysis.
  // It could also be used (cautiously) for other purposes.
  LoadGridDataAtStart = TRUE;


#ifdef USE_PYTHON
  NumberOfPythonCalls = 0;
  NumberOfPythonTopGridCalls = 0;
  NumberOfPythonSubcycleCalls = 0;
  grid_dictionary = PyDict_New();
  old_grid_dictionary = PyDict_New();
  hierarchy_information = PyDict_New();
  yt_parameter_file = PyDict_New();
  conversion_factors = PyDict_New();
  my_processor = PyLong_FromLong((Eint) MyProcessorNumber);
#endif

  /* Some stateful variables for EvolveLevel */
  for(i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++) {
    LevelCycleCount[i] = 0;
    dtThisLevelSoFar[i] = dtThisLevel[i] = 0.0;
  }


  /* Gas drag parameters */
  // UseGasDrag = 0;
  // GasDragCoefficient = 0.;

  return SUCCESS;
}
