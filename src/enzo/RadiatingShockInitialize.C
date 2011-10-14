/***********************************************************************
/
/  INITIALIZE RADIATING SEDOV BLAST WAVE
/
/  written by: Brian O'Shea
/  date:       August 2007
/  modified1:  
/
/  PURPOSE:
/
/   REFERENCE: Self-similar solution: L.I. Sedov (1946);
/              see also: Sedov (1959), Similarity and Dimensional Methods
/              in Mechanics, pp. 210, 219, 228;
/              see also: Landau & Lifshitz, Fluid Dynamics, Sect. 99
/              "The Propagation of Strong Shock Waves" (1959).
/              Experiments, terrestrial/numerical: Taylor (1941, 1949).
/
/   Two dimensional parameters: explosion energy E and ambient density rho_1
/   Two independent variables: radius r, time t
/   One dimensionless combination: r*(rho_1/E/t^2)^(1/5)
/
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
#define DEFINE_STORAGE
#include "RadiatingShockGlobalData.h"
#undef DEFINE_STORAGE

/* Set default parameter values. */

const char config_radiating_shock_defaults[] = 
"### RADIATING SHOCK DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    RadiatingShock: {\n"
"       InnerDensity  =1.0;\n"
"       OuterDensity  =1.0;\n"
"       Pressure  =1e-5;\n"
"       Energy  =1.0;\n"
"       SubgridLeft =0;\n"
"       SubgridRight  =1;\n"
"       UseDensityFluctuations  =0;\n"
"       RandomSeed  =123456789;\n"
"       DensityFluctuationLevel =0.1;\n"
"       InitializeWithKE  =0;\n"
"       SedovBlastRadius  =0.05;\n"
"       UseSedovProfile =0.05;\n"
"       KineticEnergyFraction =0.0;\n"
"       CenterPosition  =[0.5,0.5,0.5];\n"
"       SpreadOverNumZones  =3.5;\n"
"       HydrogenFractionByMass  =0.0;\n"
"       DeuteriumToHydrogenRatio  =0.0;\n"
"       HI_Fraction_Inner =0.0;\n"
"       HII_Fraction_Inner  =0.0;\n"
"       HeI_Fraction_Inner  =0.0;\n"
"       HeII_Fraction_Inner =0.0;\n"
"       HeIII_Fraction_Inner  =0.0;\n"
"       HM_Fraction_Inner =0.0;\n"
"       H2I_Fraction_Inner  =0.0;\n"
"       H2II_Fraction_Inner =0.0;\n"
"       DI_Fraction_Inner =0.0;\n"
"       DII_Fraction_Inner  =0.0;\n"
"       HDI_Fraction_Inner  =0.0;\n"
"       COI_Fraction_Inner  =0.0;\n"
"       CI_Fraction_Inner =0.0;\n"
"       CII_Fraction_Inner  =0.0;\n"
"       OI_Fraction_Inner =0.0;\n"
"       OII_Fraction_Inner  =0.0;\n"
"       SiI_Fraction_Inner  =0.0;\n"
"       SiII_Fraction_Inner =0.0;\n"
"       SiIII_Fraction_Inner  =0.0;\n"
"       CHI_Fraction_Inner  =0.0;\n"
"       CH2I_Fraction_Inner =0.0;\n"
"       CH3II_Fraction_Inner  =0.0;\n"
"       C2I_Fraction_Inner  =0.0;\n"
"       HCOII_Fraction_Inner  =0.0;\n"
"       OHI_Fraction_Inner  =0.0;\n"
"       H2OI_Fraction_Inner =0.0;\n"
"       O2I_Fraction_Inner  =0.0;\n"
"       HI_Fraction =0.0;\n"
"       HII_Fraction  =0.0;\n"
"       HeI_Fraction  =0.0;\n"
"       HeII_Fraction =0.0;\n"
"       HeIII_Fraction  =0.0;\n"
"       HM_Fraction =0.0;\n"
"       H2I_Fraction  =0.0;\n"
"       H2II_Fraction =0.0;\n"
"       DI_Fraction =0.0;\n"
"       DII_Fraction  =0.0;\n"
"       HDI_Fraction  =0.0;\n"
"       COI_Fraction  =0.0;\n"
"       CI_Fraction =0.0;\n"
"       CII_Fraction  =0.0;\n"
"       OI_Fraction =0.0;\n"
"       OII_Fraction  =0.0;\n"
"       SiI_Fraction  =0.0;\n"
"       SiII_Fraction =0.0;\n"
"       SiIII_Fraction  =0.0;\n"
"       CHI_Fraction  =0.0;\n"
"       CH2I_Fraction =0.0;\n"
"       CH3II_Fraction  =0.0;\n"
"       C2I_Fraction  =0.0;\n"
"       HCOII_Fraction  =0.0;\n"
"       OHI_Fraction  =0.0;\n"
"       H2OI_Fraction =0.0;\n"
"       O2I_Fraction  =0.0;\n"
"       UseMetallicityField =0;\n"
"       MetallicityField_Fraction =0.0;\n"
"       UseMassInjection  =0;\n"
"       InitialHydrogenMass =0.0;\n"
"       InitialDeuteriumMass  =0.0;\n"
"       InitialHeliumMass =0.0;\n"
"       InitialMetalMass  =0.0;\n"
"       MultiMetals =0;\n"
"       MultiMetalsField1_Fraction  =0.0;\n"
"       MultiMetalsField2_Fraction  =0.0;\n"
"    };\n"
"};\n";

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
 
int RadiatingShockInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
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
  char *COIName   = "COI_Density";
  char *CIName    = "CI_Density";
  char *CIIName   = "CII_Density";
  char *OIName    = "OI_Density";
  char *OIIName   = "OII_Density";
  char *SiIName   = "SiI_Density";
  char *SiIIName  = "SiII_Density";
  char *SiIIIName = "SiIII_Density";
  char *CHIName   = "CHI_Density";
  char *CH2IName  = "CH2I_Density";
  char *CH3IIName = "CH3II_Density";
  char *C2IName   = "C2I_Density";
  char *HCOIIName = "HCOII_Density";
  char *OHIName   = "OHI_Density";
  char *H2OIName  = "H2OI_Density";
  char *O2IName   = "O2I_Density";
  char *MetalName = "Metal_Density";
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

  /* parameter declarations */
 
  FLOAT RadiatingShockSubgridLeft, RadiatingShockSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT RadiatingShockCenterPosition[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, k, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 2D or 3D */
 
  if (MetaData.TopGridRank < 2 || MetaData.TopGridRank > 3) {
    ENZO_VFAIL("Cannot do RadiatingShock in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }
 
  /* There are many parameters:  geometry (cylindrical or spherical symmetry),
                                 gamma,
				 E,
				 rho_1.
     Set their default values */
 
  int RadiatingShockRandomSeedInitialize = 0;

  for(i=0; i<MAX_DIMENSION; i++)
    RadiatingShockCenterPosition[i] = 0.5;  // right in the middle of the box

  float Pi                      = 3.14159;
  float RadiatingShockVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float RadiatingShockPressure      = 1e-5;
  float RadiatingShockInnerDensity             = 1.0;
  float RadiatingShockOuterDensity             = 1.0;
  float RadiatingShockEnergy        = 1.0;
  int RadiatingShockUseDensityFluctuations = 0;
  int RadiatingShockRandomSeed = 123456789;
  float RadiatingShockDensityFluctuationLevel = 0.1;
  int RadiatingShockInitializeWithKE = 0;
  int RadiatingShockUseSedovProfile = 0;
  FLOAT RadiatingShockSedovBlastRadius = 0.05;
  float RadiatingShockKineticEnergyFraction = 0.0;

  FLOAT RadiatingShockSpreadOverNumZones = 3.5;
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
    MetaData.TopGridDims[0];
 
  /* set no subgrids by default. */
 
  RadiatingShockSubgridLeft         = 0.0;    // start of subgrid(s)
  RadiatingShockSubgridRight        = 0.0;    // end of subgrid(s)


  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)
  TestProblemData.GloverChemistryModel = GloverChemistryModel; // set this from global data (kind of a hack, but necessary)
 
  /* read input from file */

  Param.Update(config_radiating_shock_defaults);

    /* read parameters specifically for radiating shock problem*/
 
    Param.GetScalar(RadiatingShockInnerDensity,"Problem.RadiatingShock.InnerDensity");
    Param.GetScalar(RadiatingShockOuterDensity,"Problem.RadiatingShock.OuterDensity");
    Param.GetScalar(RadiatingShockPressure,"Problem.RadiatingShock.Pressure");
    Param.GetScalar(RadiatingShockEnergy,"Problem.RadiatingShock.Energy");
    Param.GetScalar(RadiatingShockSubgridLeft,"Problem.RadiatingShock.SubgridLeft");
    Param.GetScalar(RadiatingShockSubgridRight,"Problem.RadiatingShock.SubgridRight");
    Param.GetScalar(RadiatingShockUseDensityFluctuations,"Problem.RadiatingShock.UseDensityFluctuations");
    Param.GetScalar(RadiatingShockRandomSeed,"Problem.RadiatingShock.RandomSeed");
    Param.GetScalar(RadiatingShockDensityFluctuationLevel,"Problem.RadiatingShock.DensityFluctuationLevel");
    Param.GetScalar(RadiatingShockInitializeWithKE,"Problem.RadiatingShock.InitializeWithKE");
    Param.GetScalar(RadiatingShockUseSedovProfile,"Problem.RadiatingShock.UseSedovProfile");
    Param.GetScalar(RadiatingShockSedovBlastRadius,"Problem.RadiatingShock.SedovBlastRadius");

    Param.GetScalar(RadiatingShockKineticEnergyFraction,"Problem.RadiatingShock.KineticEnergyFraction");

    for (i=0; i<MAX_DIMENSION; i++){
        Param.GetArray(RadiatingShockCenterPosition[i],"Problem.RadiatingShock.CenterPosition%s",i)}

    Param.GetScalar(RadiatingShockSpreadOverNumZones,"Problem.RadiatingShock.SpreadOverNumZones");

    /* read in more general test parameters to set species, turn on color fields, etc. */
    Param.GetScalar(HydrogenFractionByMass,"Problem.RadiatingShock.HydrogenFractionByMass");
    Param.GetScalar(DeuteriumToHydrogenRatio,"Problem.RadiatingShock.DeuteriumToHydrogenRatio");

    Param.GetScalar(HI_Fraction_Inner,"Problem.RadiatingShock.HI_Fraction_Inner");
    Param.GetScalar(HII_Fraction_Inner,"Problem.RadiatingShock.HII_Fraction_Inner");
    Param.GetScalar(HeI_Fraction_Inner,"Problem.RadiatingShock.HeI_Fraction_Inner");
    Param.GetScalar(HeII_Fraction_Inner,"Problem.RadiatingShock.HeII_Fraction_Inner");
    Param.GetScalar(HeIII_Fraction_Inner,"Problem.RadiatingShock.HeIII_Fraction_Inner");
    Param.GetScalar(HM_Fraction_Inner,"Problem.RadiatingShock.HM_Fraction_Inner");
    Param.GetScalar(H2I_Fraction_Inner,"Problem.RadiatingShock.H2I_Fraction_Inner");
    Param.GetScalar(H2II_Fraction_Inner,"Problem.RadiatingShock.H2II_Fraction_Inner");

    Param.GetScalar(DI_Fraction_Inner,"Problem.RadiatingShock.DI_Fraction_Inner");
    Param.GetScalar(DII_Fraction_Inner,"Problem.RadiatingShock.DII_Fraction_Inner");
    Param.GetScalar(HDI_Fraction_Inner,"Problem.RadiatingShock.HDI_Fraction_Inner");

    Param.GetScalar(COI_Fraction_Inner,"Problem.RadiatingShock.COI_Fraction_Inner");
    Param.GetScalar(CI_Fraction_Inner,"Problem.RadiatingShock.CI_Fraction_Inner");
    Param.GetScalar(CII_Fraction_Inner,"Problem.RadiatingShock.CII_Fraction_Inner");
    Param.GetScalar(OI_Fraction_Inner,"Problem.RadiatingShock.OI_Fraction_Inner");
    Param.GetScalar(OII_Fraction_Inner,"Problem.RadiatingShock.OII_Fraction_Inner");
    Param.GetScalar(SiI_Fraction_Inner,"Problem.RadiatingShock.iI_Fraction_Inner");
    Param.GetScalar(SiII_Fraction_Inner,"Problem.RadiatingShock.SiII_Fraction_Inner");
    Param.GetScalar(SiIII_Fraction_Inner,"Problem.RadiatingShock.SiIII_Fraction_Inner");
    Param.GetScalar(CHI_Fraction_Inner,"Problem.RadiatingShock.CHI_Fraction_Inner");
    Param.GetScalar(CH2I_Fraction_Inner,"Problem.RadiatingShock.CH2I_Fraction_Inner");
    Param.GetScalar(CH3II_Fraction_Inner,"Problem.RadiatingShock.CH3II_Fraction_Inner");
    Param.GetScalar(C2I_Fraction_Inner,"Problem.RadiatingShock.C2I_Fraction_Inner");
    Param.GetScalar(HCOII_Fraction_Inner,"Problem.RadiatingShock.HCOII_Fraction_Inner");
    Param.GetScalar(OHI_Fraction_Inner,"Problem.RadiatingShock.OHI_Fraction_Inner");
    Param.GetScalar(H2OI_Fraction_Inner,"Problem.RadiatingShock.H2OI_Fraction_Inner");
    Param.GetScalar(O2I_Fraction_Inner,"Problem.RadiatingShock.O2I_Fraction_Inner");

    Param.GetScalar(HI_Fraction,"Problem.RadiatingShock.HI_Fraction");
    Param.GetScalar(HII_Fraction,"Problem.RadiatingShock.HII_Fraction");
    Param.GetScalar(HeI_Fraction,"Problem.RadiatingShock.HeI_Fraction");
    Param.GetScalar(HeII_Fraction,"Problem.RadiatingShock.HeII_Fraction");
    Param.GetScalar(HeIII_Fraction,"Problem.RadiatingShock.HeIII_Fraction");
    Param.GetScalar(HM_Fraction,"Problem.RadiatingShock.HM_Fraction");
    Param.GetScalar(H2I_Fraction,"Problem.RadiatingShock.H2I_Fraction");
    Param.GetScalar(H2II_Fraction,"Problem.RadiatingShock.H2II_Fraction");

    Param.GetScalar(DI_Fraction,"Problem.RadiatingShock.DI_Fraction");
    Param.GetScalar(DII_Fraction,"Problem.RadiatingShock.DII_Fraction");
    Param.GetScalar(HDI_Fraction,"Problem.RadiatingShock.HDI_Fraction");

    Param.GetScalar(COI_Fraction,"Problem.RadiatingShock.COI_Fraction");
    Param.GetScalar(CI_Fraction,"Problem.RadiatingShock.CI_Fraction");
    Param.GetScalar(CII_Fraction,"Problem.RadiatingShock.CII_Fraction");
    Param.GetScalar(OI_Fraction,"Problem.RadiatingShock.OI_Fraction");
    Param.GetScalar(OII_Fraction,"Problem.RadiatingShock.OII_Fraction");
    Param.GetScalar(SiI_Fraction,"Problem.RadiatingShock.SiI_Fraction");
    Param.GetScalar(SiII_Fraction,"Problem.RadiatingShock.SiII_Fraction");
    Param.GetScalar(SiIII_Fraction,"Problem.RadiatingShock.SiIII_Fraction");
    Param.GetScalar(CHI_Fraction,"Problem.RadiatingShock.CHI_Fraction");
    Param.GetScalar(CH2I_Fraction,"Problem.RadiatingShock.CH2I_Fraction");
    Param.GetScalar(CH3II_Fraction,"Problem.RadiatingShock.CH3II_Fraction");
    Param.GetScalar(C2I_Fraction,"Problem.RadiatingShock.C2I_Fraction");
    Param.GetScalar(HCOII_Fraction,"Problem.RadiatingShock.HCOII_Fraction");
    Param.GetScalar(OHI_Fraction,"Problem.RadiatingShock.OHI_Fraction");
    Param.GetScalar(H2OI_Fraction,"Problem.RadiatingShock.H2OI_Fraction");
    Param.GetScalar(O2I_Fraction,"Problem.RadiatingShock.O2I_Fraction");

    Param.GetScalar(UseMetallicityField,"Problem.RadiatingShock.UseMetallicityField");
    Param.GetScalar(MetallicityField_Fraction,"Problem.RadiatingShock.MetallicityField_Fraction");

    Param.GetScalar(UseMassInjection,"Problem.RadiatingShock.UseMassInjection");
    Param.GetScalar(InitialHydrogenMass,"Problem.RadiatingShock.InitialHydrogenMass");
    Param.GetScalar(InitialDeuteriumMass,"Problem.RadiatingShock.InitialDeuteriumMass");
    Param.GetScalar(InitialHeliumMass,"Problem.RadiatingShock.InitialHeliumMass");
    Param.GetScalar(InitialMetalMass,"Problem.RadiatingShock.InitialMetalMass");

    Param.GetScalar(MultiMetals,"Problem.RadiatingShock.MultiMetals");
    Param.GetScalar(MultiMetalsField1_Fraction,"Problem.RadiatingShock.MultiMetalsField1_Fraction");
    Param.GetScalar(MultiMetalsField2_Fraction,"Problem.RadiatingShock.MultiMetalsField2_Fraction");

 
  /* set number of zones on the finest level to resolve the initial explosion */
  FLOAT dr = RadiatingShockSpreadOverNumZones*dx*POW(RefineBy,-MaximumRefinementLevel);
 
  /* compute p_2 as a function of explosion energy E, initial explosion
     radius dr, and gamma.
 
     2D:  p_2 = (gamma-1)*E/(pi*r^2)*rho_in
     3D:  p_2 = (gamma-1)*E/(4/3*pi*r^3)*rho_in
          rho_2 = 1 = rho_1
 
     If p_2 (inner) is too high, the code crashes in euler.src (dnu<0) after
     a few time steps.
     In 2D it is stable for RadiatingShockEnergy <= 0.025 (tested with uniform grid 200^2).
  */
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  double RadiatingShockEnergyDouble;


  // calculate some values if we aren't using the Sedov profile.  If we are, this is all
  // calculated inside Grid::RadiatingShockInitialize
  if(!RadiatingShockUseSedovProfile){
    // input is in units of 1 FOE, this converts to CGS
    RadiatingShockEnergyDouble = double(RadiatingShockEnergy) * 1.0e+51;

    // converts this to ergs per gram in code units
    RadiatingShockEnergyDouble /= (double(DensityUnits) * POW(double(LengthUnits),5.0) * POW(double(TimeUnits),-2.0) );
  
    RadiatingShockEnergy = float(RadiatingShockEnergyDouble);

  } else {
    printf("using sedov profile: using RadiatingShockEnergy alone!\n");
 
  }

  if(debug)
    printf("RadiatingShockInitialize:  RadiatingShockEnergy is %e in code units\n",RadiatingShockEnergy);

  float numberInjectionCells = 0;
  double InjectionMass2Density_scaleFactor;

  // If injecting a mass of gas into the center, 
  // effective number of cells is simply area/volume of circle/sphere
  // with r = RadiatingShockSpreadOverNumZones.
  if (TestProblemData.UseMassInjection) {
    // Make sure total mass is not zero.
    if ((TestProblemData.InitialHydrogenMass <= 0.0) && (TestProblemData.InitialHeliumMass <= 0.0)) {
      ENZO_FAIL("Hydrogen and helium mass cannot both be zero.  That would be zero mass in the center.\n");
    }

    // 2D
    if(MetaData.TopGridRank==2){
      numberInjectionCells = Pi * pow(RadiatingShockSpreadOverNumZones,2);
    }
    // 3D
    else {
      numberInjectionCells = (4./3.) * Pi * pow(RadiatingShockSpreadOverNumZones,3);
    }

    InjectionMass2Density_scaleFactor = 1.989e33 / 
      (MassUnits * numberInjectionCells * 
       pow((dx*POW(RefineBy,-MaximumRefinementLevel)),3));

    // ignore D mass
    RadiatingShockInnerDensity = (TestProblemData.InitialHydrogenMass + 
				  TestProblemData.InitialHeliumMass) * 
      float(InjectionMass2Density_scaleFactor);
    fprintf(stderr,"Setting inner density to %e.\n",RadiatingShockInnerDensity);

    // Adjust H mass fraction
    TestProblemData.InnerHydrogenFractionByMass = TestProblemData.InitialHydrogenMass /
      (TestProblemData.InitialHydrogenMass + TestProblemData.InitialHeliumMass);
    fprintf(stderr,"Inner H mass fraction is %.2f.\n",TestProblemData.InnerHydrogenFractionByMass);

    // Adjust D/H ratio
    TestProblemData.InnerDeuteriumToHydrogenRatio = TestProblemData.InitialDeuteriumMass /
      TestProblemData.InitialHydrogenMass;
    fprintf(stderr,"Inner D/H ratio is %.2e.\n",TestProblemData.InnerDeuteriumToHydrogenRatio);

    // Set metal fraction
    TestProblemData.MetallicityField_Fraction = TestProblemData.InitialMetalMass /
      (TestProblemData.InitialHydrogenMass + TestProblemData.InitialHeliumMass);
    fprintf(stderr,"Inner metal fraction is %.2e.\n",TestProblemData.MetallicityField_Fraction);
  }
  else {
    TestProblemData.InnerHydrogenFractionByMass = TestProblemData.HydrogenFractionByMass;
    TestProblemData.InnerDeuteriumToHydrogenRatio = TestProblemData.DeuteriumToHydrogenRatio;
  }

  // BWO: modified to include the idea that the inner region might have more gas
  float RadiatingShockInnerPressure = 3.0*(Gamma-1.0)*RadiatingShockEnergy*RadiatingShockInnerDensity/
                                  (MetaData.TopGridRank + 1.0)/
                                  POW(dr,MetaData.TopGridRank)/Pi;
 
  /* Check the self-similarity condition: p2/p1 >> (gamma+1)/(gamma-1). */
 
  float pjump = RadiatingShockInnerPressure/RadiatingShockPressure;
  if ( pjump < 10.*(Gamma+1)/(Gamma-1) )
    printf("SBI: WARNING! No self-similarity. Pressure jump %"FSYM".\n", pjump);
 
  /* Compute energies */

  // ambient gas internal energy
  RadiatingShockTotalEnergy = RadiatingShockPressure/((Gamma - 1.0)*RadiatingShockOuterDensity);

  if(debug)
    printf("ambient gas energy should be %e\n",RadiatingShockTotalEnergy);

  // heated gas internal energy
  RadiatingShockInnerTotalEnergy= RadiatingShockInnerPressure/((Gamma - 1.0)*
						   RadiatingShockInnerDensity);

  /* calculate kinetic energy quantities (if we're initializing with a sawtooth) 
     this is similar to what we do above, for 2D (cylindrical) and 3D(spherical) */
  float RadiatingShockRhoZero, RadiatingShockVelocityZero,  MassZero;

  RadiatingShockRhoZero = RadiatingShockVelocityZero =  MassZero = 0.0;

  /* we calculate all of this stuff when not using a Sedov blast profile.  If we ARE using
     Sedov values, ignore it because it doesn't matter. */
  if(RadiatingShockInitializeWithKE && !RadiatingShockUseSedovProfile){

    if(MetaData.TopGridRank==2){  // cylindrical supernova (2D)

      MassZero = Pi * POW(dr,2.0) * RadiatingShockInnerDensity; // Mzero in code units (actually mass per unit length)

      RadiatingShockRhoZero = MassZero / (2.0*Pi/3.0* POW(dr,2.0) );  // this is really density

      // units actually work out to velocity (assuming input energy is really energy per unit length)
      RadiatingShockVelocityZero = POW( 10.0/3.0 * (RadiatingShockEnergy*RadiatingShockKineticEnergyFraction) / MassZero, 0.5 );

    } else {  // spherical supernova (3D)

      MassZero = 1.3333 * Pi * POW(dr,3.0) * RadiatingShockInnerDensity;  // Mzero in code units

      RadiatingShockRhoZero = MassZero / Pi / POW(dr,3.0);

      fprintf(stderr,"SBI: RadiatingShockKineticEnergyFraction is %e\n",RadiatingShockKineticEnergyFraction);

      RadiatingShockVelocityZero = POW( 3.0 * (RadiatingShockEnergy*RadiatingShockKineticEnergyFraction) / MassZero, 0.5);

    }

  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }
 
  /* set up uniform top grid before setting up explosion */
 
  if (TopGrid.GridData->InitializeUniformGrid(RadiatingShockOuterDensity,
					      RadiatingShockTotalEnergy,
					      RadiatingShockTotalEnergy,
					      RadiatingShockVelocity,
                          ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
  /* Create as many subgrids as there are refinement levels
     needed to resolve the initial explosion region upon the start-up. */
 
  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0)
    Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
 
  /* Create new HierarchyEntries. */
 
  int lev;
  for (lev = 0; lev < MaximumRefinementLevel; lev++)
    Subgrid[lev] = new HierarchyEntry;

  for (lev = 0; lev < MaximumRefinementLevel; lev++) {

    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((RadiatingShockSubgridRight - RadiatingShockSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("RadiatingShock:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
	     NumberOfSubgridZones[0]);
 
    if (NumberOfSubgridZones[0] > 0) {
 
      /* fill them out */
 
      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == MaximumRefinementLevel-1)
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
	LeftEdge[dim]    = RadiatingShockSubgridLeft;
	RightEdge[dim]   = RadiatingShockSubgridRight;
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(RadiatingShockOuterDensity,
						   RadiatingShockTotalEnergy,
						   RadiatingShockTotalEnergy,
						   RadiatingShockVelocity,
                           ZeroBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->RadiatingShockInitializeGrid(dr,
				    RadiatingShockInnerDensity,
				    RadiatingShockInnerTotalEnergy,
				    RadiatingShockUseDensityFluctuations,
				    RadiatingShockRandomSeed,
				    RadiatingShockDensityFluctuationLevel,
				    RadiatingShockInitializeWithKE,
				    RadiatingShockUseSedovProfile,
				    RadiatingShockSedovBlastRadius,
				    RadiatingShockEnergy,
				    RadiatingShockPressure,
				    RadiatingShockKineticEnergyFraction,
				    RadiatingShockRhoZero,
				    RadiatingShockVelocityZero,
				    RadiatingShockRandomSeedInitialize,
 				    RadiatingShockCenterPosition) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in RadiatingShockInitialize[Sub]Grid.");
	}

      RadiatingShockRandomSeedInitialize = 1;  // random number generator is now seeded - don't do it for topgrid
    }
    else{
      printf("RadiatingShock: single grid start-up.\n");
    }
  }

  /* set up subgrids from level 1 to max refinement level -1 */
 
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  /* set up the root grid */
 
  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (TopGrid.GridData->RadiatingShockInitializeGrid(dr,
				    RadiatingShockInnerDensity,
			            RadiatingShockInnerTotalEnergy,
				    RadiatingShockUseDensityFluctuations,
				    RadiatingShockRandomSeed,
				    RadiatingShockDensityFluctuationLevel,
				    RadiatingShockInitializeWithKE,
				    RadiatingShockUseSedovProfile,
				    RadiatingShockSedovBlastRadius,
				    RadiatingShockEnergy,
				    RadiatingShockPressure,
				    RadiatingShockKineticEnergyFraction,
				    RadiatingShockRhoZero,
				    RadiatingShockVelocityZero,
				    RadiatingShockRandomSeedInitialize,
 				    RadiatingShockCenterPosition ) == FAIL) {
            ENZO_FAIL("Error in RadiatingShockInitializeGrid.");
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

  if(TestProblemData.GloverChemistryModel){

    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

    DataLabel[i++] = HIIName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = H2IName;

    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
      DataLabel[i++] = HeIName;
      DataLabel[i++] = HeIIName;
      DataLabel[i++] = HeIIIName;
    }

    if( (GCM==3) || (GCM==5) || (GCM==7) ){
      DataLabel[i++] = COIName;
    }

    if( (GCM==2) || (GCM==3) || (GCM==7) ){
      DataLabel[i++] = CIName;
      DataLabel[i++] = CIIName;
      DataLabel[i++] = OIName;
      DataLabel[i++] = OIIName;
    }

    if( (GCM==2) || (GCM==3) ){
      DataLabel[i++] = SiIName;
      DataLabel[i++] = SiIIName;
      DataLabel[i++] = SiIIIName;
    }

    if( (GCM==3) || (GCM==7) ){
      DataLabel[i++] = CHIName;
      DataLabel[i++] = CH2IName;
      DataLabel[i++] = CH3IIName;
      DataLabel[i++] = C2IName;
      DataLabel[i++] = HCOIIName;
      DataLabel[i++] = OHIName;
      DataLabel[i++] = H2OIName;
      DataLabel[i++] = O2IName;
    }

  } //   if(TestProblemData.GloverChemistryModel)


  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "RadiatingShockInnerDensity         = %"FSYM"\n"  , RadiatingShockInnerDensity);
    fprintf(Outfptr, "RadiatingShockOuterDensity         = %"FSYM"\n"  , RadiatingShockOuterDensity);
    fprintf(Outfptr, "RadiatingShockPressure        = %"FSYM"\n"  , RadiatingShockPressure);
    fprintf(Outfptr, "RadiatingShockEnergy          = %"FSYM"\n"  , RadiatingShockEnergy);
    fprintf(Outfptr, "RadiatingShockInnerPressure   = %"FSYM"\n"  ,
	    RadiatingShockInnerPressure);
 
    fprintf(Outfptr, "RadiatingShockUseDensityFluctuations   = %"ISYM"\n", RadiatingShockUseDensityFluctuations);
    fprintf(Outfptr, "RadiatingShockRandomSeed   = %"ISYM"\n", RadiatingShockRandomSeed);
    fprintf(Outfptr, "RadiatingShockDensityFluctuationLevel   = %"FSYM"\n", RadiatingShockDensityFluctuationLevel);
    fprintf(Outfptr, "RadiatingShockInitializeWithKE = %"ISYM"\n", RadiatingShockInitializeWithKE);
    fprintf(Outfptr, "RadiatingShockUseSedovProfile = %"ISYM"\n", RadiatingShockUseSedovProfile);

    fprintf(Outfptr, "RadiatingShockSedovBlastRadius = %"PSYM"\n", RadiatingShockSedovBlastRadius);

    fprintf(Outfptr,  "RadiatingShockKineticEnergyFraction = %"FSYM"\n", RadiatingShockKineticEnergyFraction);

    fprintf(Outfptr, "RadiatingShockCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RadiatingShockCenterPosition, RadiatingShockCenterPosition+1,
		  RadiatingShockCenterPosition+2);

    fprintf(Outfptr, "RadiatingShockSpreadOverNumZones  = %"PSYM"\n", RadiatingShockSpreadOverNumZones);

    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);

    fprintf(Outfptr, "TestProblemInitialHIFractionInner  = %"FSYM"\n", TestProblemData.HI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHIIFractionInner  = %"FSYM"\n", TestProblemData.HII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIFractionInner  = %"FSYM"\n", TestProblemData.HeI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIFractionInner  = %"FSYM"\n", TestProblemData.HeII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIIFractionInner  = %"FSYM"\n", TestProblemData.HeIII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHMFractionInner  = %"FSYM"\n", TestProblemData.HM_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IFractionInner  = %"FSYM"\n", TestProblemData.H2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IIFractionInner  = %"FSYM"\n", TestProblemData.H2II_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialDIFractionInner  = %"FSYM"\n", TestProblemData.DI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialDIIFractionInner  = %"FSYM"\n", TestProblemData.DII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHDIFractionInner  = %"FSYM"\n", TestProblemData.HDI_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialCOIFractionInner  = %"FSYM"\n", TestProblemData.COI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCIFractionInner  = %"FSYM"\n", TestProblemData.CI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCIIFractionInner  = %"FSYM"\n", TestProblemData.CII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOIFractionInner  = %"FSYM"\n", TestProblemData.OI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOIIFractionInner  = %"FSYM"\n", TestProblemData.OII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIFractionInner  = %"FSYM"\n", TestProblemData.SiI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIIFractionInner  = %"FSYM"\n", TestProblemData.SiII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIIIFractionInner  = %"FSYM"\n", TestProblemData.SiIII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCHIFractionInner  = %"FSYM"\n", TestProblemData.CHI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCH2IFractionInner  = %"FSYM"\n", TestProblemData.CH2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCH3IIFractionInner  = %"FSYM"\n", TestProblemData.CH3II_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialC2IFractionInner  = %"FSYM"\n", TestProblemData.C2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHCOIIFractionInner  = %"FSYM"\n", TestProblemData.HCOII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOHIFractionInner  = %"FSYM"\n", TestProblemData.OHI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2OIFractionInner  = %"FSYM"\n", TestProblemData.H2OI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialO2IFractionInner  = %"FSYM"\n", TestProblemData.O2I_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);

    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);

    fprintf(Outfptr, "TestProblemInitialCOIFraction  = %"FSYM"\n", TestProblemData.COI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCIFraction  = %"FSYM"\n", TestProblemData.CI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCIIFraction  = %"FSYM"\n", TestProblemData.CII_Fraction);
    fprintf(Outfptr, "TestProblemInitialOIFraction  = %"FSYM"\n", TestProblemData.OI_Fraction);
    fprintf(Outfptr, "TestProblemInitialOIIFraction  = %"FSYM"\n", TestProblemData.OII_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIFraction  = %"FSYM"\n", TestProblemData.SiI_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIIFraction  = %"FSYM"\n", TestProblemData.SiII_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIIIFraction  = %"FSYM"\n", TestProblemData.SiIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialCHIFraction  = %"FSYM"\n", TestProblemData.CHI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCH2IFraction  = %"FSYM"\n", TestProblemData.CH2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialCH3IIFraction  = %"FSYM"\n", TestProblemData.CH3II_Fraction);
    fprintf(Outfptr, "TestProblemInitialC2IFraction  = %"FSYM"\n", TestProblemData.C2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialHCOIIFraction  = %"FSYM"\n", TestProblemData.HCOII_Fraction);
    fprintf(Outfptr, "TestProblemInitialOHIFraction  = %"FSYM"\n", TestProblemData.OHI_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2OIFraction  = %"FSYM"\n", TestProblemData.H2OI_Fraction);
    fprintf(Outfptr, "TestProblemInitialO2IFraction  = %"FSYM"\n", TestProblemData.O2I_Fraction);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

    fprintf(Outfptr, "TestProblemUseMassInjection  = %"ISYM"\n", TestProblemData.UseMassInjection);
    fprintf(Outfptr, "TestProblemInitialHydrogenMass  = %"ESYM"\n", TestProblemData.InitialHydrogenMass);
    fprintf(Outfptr, "TestProblemInitialDeuteriumMass  = %"ESYM"\n", TestProblemData.InitialDeuteriumMass);
    fprintf(Outfptr, "TestProblemInitialHeliumMass  = %"ESYM"\n", TestProblemData.InitialHeliumMass);
    fprintf(Outfptr, "TestProblemInitialMetalMass  = %"ESYM"\n", TestProblemData.InitialMetalMass);

    fprintf(Outfptr, "TestProblemMultiMetals  = %"ISYM"\n", TestProblemData.MultiMetals);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField1_Fraction);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField2_Fraction);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 

 
  return SUCCESS;
 
}
