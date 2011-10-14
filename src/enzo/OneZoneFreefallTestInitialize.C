/***********************************************************************
/
/  INITIALIZE ONE-ZONE FREE-FALL TEST PROBLEM
/
/  written by: Britton Smith
/  date:       JANUARY 2009
/  modified1:  
/
/  PURPOSE: Test chemistry and cooling in free-fall collapse.
/
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

/* Set default parameter values. */

const char config_one_zone_freefall_test_defaults[] =
"### ONE ZONE FREEFALL TEST DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    OneZoneFreefallTest: {\n"
"        InitialDensity 	= 1.0;\n"
"        MinimumEnergy 		= 10.0;\n"
"        MaximumEnergy 		= 1000.0;\n"
"        MinimumMetallicity 	= 1e-6;\n"
"        MaximumMetallicity 	= 1e-2;\n"
"        ConstantDensityVelocity	= [0.0,0.0,0.0];\n"
"        MaximumRefinementLevel = 0;\n"
"\n"
"        InitialHIFraction 	= -99999.9;\n"
"        InitialHIIFraction	= -99999.9;\n" 
"        InitialHeIFraction	= -99999.9;\n"
"        InitialHeIIFraction	= -99999.9;\n"
"        InitialHeIIIFraction 	= -99999.9;\n"
"        InitialHMFraction 	= -99999.9;\n"
"        InitialH2IFraction	= -99999.9;\n"
"        InitialH2IIFraction 	= -99999.9;\n"
"        InitialDIFraction 	= -99999.9;\n"
"        InitialDIIFraction 	= -99999.9;\n"
"        InitialHDIFraction 	= -99999.9;\n"
"        UseMetallicityField 	= FALSE;\n"
"        MetallicityNormalization  = 1.0\n"
"        \n"
"    };\n"
"};\n";


int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int OneZoneFreefallTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
				  TopGridData &MetaData)
{

  fprintf(stderr,"Initializing one-zone free-fall test.\n");

  char *DensName = "Density";
  char *TEName   = "Total_Energy";
  char *GEName   = "Gas_Energy";
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
  char *MetalName = "Metal_Density";
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

   /* parameter declarations */
 
  FLOAT ConstantDensitySubgridLeft, ConstantDensitySubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  float ConstantDensityVelocity[3];

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

  MaximumRefinementLevel = 0;

  float OneZoneFreefallTestInitialDensity;
  float OneZoneFreefallTestMinimumEnergy;
  float OneZoneFreefallTestMaximumEnergy;
  float OneZoneFreefallTestMinimumMetallicity;
  float OneZoneFreefallTestMaximumMetallicity;
  TestProblemData.OneZoneFreefallTimestepFraction;

//   /* set no subgrids by default. */
 
  ConstantDensitySubgridLeft         = 0.0;    // start of subgrid(s)
  ConstantDensitySubgridRight        = 0.0;    // end of subgrid(s)

  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_one_zone_freefall_test_defaults);


  /* read parameters */

  /* read in more general test parameters to set species, turn on color fields, etc. */
  Param.GetScalar(OneZoneFreefallTestInitialDensity,
		"Problem.OneZoneFreefallTest.InitialDensity");
  Param.GetScalar(OneZoneFreefallTestMinimumEnergy,
		"Problem.OneZoneFreefallTest.MinimumEnergy");
  Param.GetScalar(OneZoneFreefallTestMaximumEnergy,
		"Problem.OneZoneFreefallTest.MaximumEnergy");
  Param.GetScalar(OneZoneFreefallTestMinimumMetallicity,
		"Problem.OneZoneFreefallTest.MinimumMetallicity");
  Param.GetScalar(OneZoneFreefallTestMaximumMetallicity,
		"Problem.OneZoneFreefallTest.MaximumMetallicity");
  Param.GetScalar(OneZoneFreefallTimestepFraction,
		"Problem.OneZoneFreefall.TimestepFraction");

  Param.GetScalar(TestProblemHydrogenFractionByMass,
		"Problem.OneZoneFreefall.HydrogenFractionByMass");
  Param.GetScalar(TestProblemDeuteriumToHydrogenRatio,
		"Problem.OneZoneFreefall.DeuteriumToHydrogenRatio");
  Param.GetScalar(TestProblemInitialHIFraction,
		"Problem.OneZoneFreefall.HI_Fraction");
  Param.GetScalar(TestProblemInitialHIIFraction, 
		"Problem.OneZoneFreefall.HII_Fraction");
  Param.GetScalar(TestProblemInitialHeIFraction, 
		"Problem.OneZoneFreefall.HeI_Fraction");
  Param.GetScalar(TestProblemInitialHeIIFraction, 
		"Problem.OneZoneFreefall.HeII_Fraction");
  Param.GetScalar(TestProblemInitialHeIIIIFraction, 
		"Problem.OneZoneFreefall.HeIII_Fraction");
  Param.GetScalar(TestProblemInitialHMFraction, 
		"Problem.OneZoneFreefall.HM_Fraction");
  Param.GetScalar(TestProblemInitialH2IFraction, 
		"Problem.OneZoneFreefall.H2I_Fraction");
  Param.GetScalar(TestProblemInitialH2IIFraction, 
		"Problem.OneZoneFreefall.H2II_Fraction");
  Param.GetScalar(TestProblemInitialDIFraction, 
		"Problem.OneZoneFreefall.DI_Fraction");
  Param.GetScalar(TestProblemInitialDIIFraction, 
		"Problem.OneZoneFreefall.DII_Fraction");
  Param.GetScalar(TestProblemInitialHDIFraction, 
		"Problem.OneZoneFreefall.HDI_Fraction");
  Param.GetScalar(TestProblemUseMetallicityField,
		"Problem.OneZoneFreefall.UseMetallicityField");
  Param.GetScalar(TestProblemMetallicityNormalization, 
		"Problem.OneZoneFreefall.MetallicityNormalization");

  /* Set constant for analytical free-fall collapse. */
  TestProblemData.OneZoneFreefallConstant = pow(OneZoneFreefallTestInitialDensity, -0.5);

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, 0.0) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }
 
 //  /* set up grid */
 
  if (TopGrid.GridData->OneZoneFreefallTestInitializeGrid(OneZoneFreefallTestInitialDensity,
							  OneZoneFreefallTestMinimumEnergy,
							  OneZoneFreefallTestMaximumEnergy,
							  OneZoneFreefallTestMinimumMetallicity,
							  OneZoneFreefallTestMaximumMetallicity) == FAIL) {
    ENZO_FAIL("Error in OneZoneFreefallTestInitializeGrid.");
  }
 
  /* set up field names and units -- NOTE: these absolutely MUST be in 
     the same order that they are in Grid_InitializeUniformGrids.C, or 
     else you'll find out that data gets written into incorrectly-named
     fields.  Just FYI. */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;

  if(MetaData.TopGridRank > 1)
    DataLabel[i++] = Vel2Name;

  if(MetaData.TopGridRank > 2)
    DataLabel[i++] = Vel3Name;

  if (TestProblemData.MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }
 
  if (TestProblemData.UseMetallicityField) {
    DataLabel[i++] = MetalName;

    if(TestProblemData.MultiMetals){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
  }
 
  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "OneZoneFreefallTestInitialDensity = %"FSYM, OneZoneFreefallTestInitialDensity);
    fprintf(Outfptr, "OneZoneFreefallTestMinimumEnergy = %"FSYM, OneZoneFreefallTestMinimumEnergy);
    fprintf(Outfptr, "OneZoneFreefallTestMaximumEnergy = %"FSYM, OneZoneFreefallTestMaximumEnergy);
    fprintf(Outfptr, "OneZoneFreefallTestMinimumMetallicity = %"FSYM, OneZoneFreefallTestMinimumMetallicity);
    fprintf(Outfptr, "OneZoneFreefallTestMaximumMetallicity = %"FSYM, OneZoneFreefallTestMaximumMetallicity);
    fprintf(Outfptr, "OneZoneFreefallTimestepFraction = %"FSYM, TestProblemData.OneZoneFreefallTimestepFraction);

    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);
    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemMetallicityNormalization  = %"FSYM"\n", TestProblemData.MetallicityNormalization);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 
 
  return SUCCESS;
 
}
