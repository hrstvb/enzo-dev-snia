/***********************************************************************
/
/  INITIALIZE COOLING TEST PROBLEM
/
/  written by: Britton Smith
/  date:       JANUARY 2009
/  modified1:  
/
/  PURPOSE: Test cooling methods without hydrodynamics.
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
   
    
const char config_cooling_test_defaults[] = 
"### COOOLING TEST DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    Cooling Test: {\n"
"       MinimumHNumberDensity    = 1.0;\n"   
"       MaximumHNumberDensity    = 1e6;\n"  
"       MinimumMetallicity       = 1e-6;\n"	     
"       MaximumMetallicity       = 1;\n"	     
"       MinimumTemperature       = 10;\n"	     
"       MaximumTemperature       = 1e7;\n"	     
"       ResetEnergies            = 1;\n"
"\n"	     
"       HydrogenFractionByMass   = 0.76;\n"  
"       DeuteriumToHydrogenRatio = 2.0*3.4e-5 ; # Burles & Tytler 1998 \n"
"\n"
"       HI_Fraction              = 1.0;\n"	     
"       HII_Fraction             = tiny_number;\n"	     
"       HeI_Fraction             = 1.0;\n"	     
"       HeIIFraction             = tiny_number;\n"	     
"       HeIII_Fraction           = tiny_number;\n"	     
"       HM_Fraction              = tiny_number;\n"	     
"       H2I_Fraction             = tiny_number;\n"	     
"       H2II_Fraction            = tiny_number;\n"	     
"       UseMetallicityField      = False;\n"     
"\n"  
"    };\n"
"};\n";

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int CoolingTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData)
{

  fprintf(stderr,"Initializing cooling test.\n");

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

  float ConstantDensityVelocity[3]   = {0.0, 0.0, 0.0};

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

  MaximumRefinementLevel = 0;

  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
     MetaData.TopGridDims[0];
 
//   /* set no subgrids by default. */
 
  ConstantDensitySubgridLeft         = 0.0;    // start of subgrid(s)
  ConstantDensitySubgridRight        = 0.0;    // end of subgrid(s)

  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

  /* read input from file */

  int comment_count = 0;
 
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {
 
    ret = 0;
 
    /* read parameters specifically for constant density problem */

    /* read in more general test parameters to set species, turn on color fields, etc. */
    Param.GetScalar(TestProblemData.MinimumHNumberDensity, "Problem.CoolingTest.MinimumHNumberDensity");
    Param.GetScalar(TestProblemData.MaximumHNumberDensity,"Problem.CoolingTest.MaximumHNumberDensity ");
    Param.GetScalar(TestProblemData.MinimumMetallicity,"Problem.CoolingTest.MinimumMetallicity");
    Param.GetScalar(TestProblemData.MaximumMetallicity,"Problem.CoolingTest.MaximumMetallicity");
    Param.GetScalar(TestProblemData.MinimumTemperature,"Problem.CoolingTest.MinimumTemperature");
    Param.GetScalar(TestProblemData.MaximumTemperature,"Problem.CoolingTest.MaximumTemperature");
    Param.GetScalar(TestProblemData.ResetEnergies,"Problem.CoolingTest.ResetEnergies");
    Param.GetScalar(TestProblemData.HydrogenFractionByMass,"Problem.CoolingTest.HydrogenFractionByMass");
    Param.GetScalar(TestProblemData.DeuteriumToHydrogenRatio,"Problem.CoolingTest.DeuteriumToHydrogenRatio");
    Param.GetScalar(TestProblemData.HI_Fraction,"Problem.CoolingTest.HI_Fraction ");
    Param.GetScalar(TestProblemData.HII_Fraction,"Problem.CoolingTest.HII_Fraction");
    Param.GetScalar(TestProblemData.HeI_Fraction,"Problem.CoolingTest.HeI_Fraction");
    Param.GetScalar(TestProblemData.HeII_Fraction,"Problem.CoolingTest.HeIIFraction");
    Param.GetScalar(TestProblemData.HeIII_Fraction,"Problem.CoolingTest.HeIII_Fraction");
    Param.GetScalar(TestProblemData.HM_Fraction,"Problem.CoolingTest.HM_Fraction");
    Param.GetScalar(TestProblemData.H2I_Fraction,"Problem.CoolingTest.H2I_Fraction");
    Param.GetScalar(TestProblemData.H2II_Fraction,"Problem.CoolingTest.H2II_Fraction");		   
    Param.GetScalar(TestProblemData.UseMetallicityField,"Problem.CoolingTest.UseMetallicityField");
    

 
  } // end input from parameter file

  // Use metallicity field.
  TestProblemData.UseMetallicityField = 1;

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, 0.0) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }
 
 //  /* set up grid */
 
  if (TopGrid.GridData->CoolingTestInitializeGrid() == FAIL) {
    fprintf(stderr, "Error in InitializeCoolingTestGrid.\n");
    return FAIL;
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
    fprintf(Outfptr, "CoolingTestMinimumHNumberDensity = %"FSYM, TestProblemData.MinimumHNumberDensity);
    fprintf(Outfptr, "CoolingTestMaximumHNumberDensity = %"FSYM, TestProblemData.MaximumHNumberDensity);
    fprintf(Outfptr, "CoolingTestMinimumMetallicity = %"FSYM, TestProblemData.MinimumMetallicity);
    fprintf(Outfptr, "CoolingTestMaximumMetallicity = %"FSYM, TestProblemData.MaximumMetallicity);
    fprintf(Outfptr, "CoolingTestMinimumTemperature = %"FSYM, TestProblemData.MinimumTemperature);
    fprintf(Outfptr, "CoolingTestMaximumTemperature = %"FSYM, TestProblemData.MaximumTemperature);

    fprintf(Outfptr, "CoolingTestResetEnergies = %"ISYM, TestProblemData.ResetEnergies);

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

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 
 
  return SUCCESS;
 
}
