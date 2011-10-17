/***********************************************************************
/
/  GRID CLASS (COMPUTE ELECTRON FRACTION ESTIMATE)
/
/  written by: John H. Wise
/  date:       April, 2007
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#define EFRAC_LOWERLIMIT 1e-2
#define MIN_TEMP 50
#define MAX_IONTEMP 30000

#include <stdlib.h>
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
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::ElectronFractionEstimate(float dt)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Find necessary fields. */

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, gammaNum;

  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
			     Vel3Num, TENum);

  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, 
        VelocityUnits, TimeUnits, aUnits = 1;
  FLOAT a = 1.0, dadt;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  float afloat = float(a);

  /* For cells with photo-ionization rates and low e-fractions (shell
     of HII regions), estimate e-fraction. */

  float mh = 1.673e-24;
  float alpha_recombination = 2.59e-13;  // Assume T = 1e4 K
  alpha_recombination *= (TimeUnits * (DensityUnits / mh));

  int i, j, k, index;
  float efrac, t_i, x_eq, x_estimate;
  float total, total_h, t_frac, new_hii;

  // Compton cooling
  float comp1 = 1e-20, comp2 = 1e-20, zr;
  if (ComovingCoordinates) {
    zr = 1.0 / (afloat * aUnits) - 1.0;
    comp1 = CoolData.comp * (1.0 + zr);
    comp2 = 2.73 * (1.0 + zr);
  }

  double CoolUnit, xbase1, dbase1;
  xbase1 = LengthUnits / (afloat * aUnits);
  dbase1 = DensityUnits * pow((afloat*aUnits), 3);
  CoolUnit = (pow(aUnits,5) * pow(xbase1,2) * pow(mh,2)) /
    (pow(TimeUnits,3) * dbase1);

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	efrac = BaryonField[DeNum][index] / BaryonField[DensNum][index];

	if (BaryonField[kphHINum][index] > 0 && efrac < EFRAC_LOWERLIMIT) {

	  /* Estimate ionization fraction */

	  total = BaryonField[kphHINum][index] + 
	    (alpha_recombination * BaryonField[DeNum][index]);
	  t_i = 1.0 / total;
	  x_eq = BaryonField[kphHINum][index] / total;
	  t_frac = dt / t_i;
	  x_estimate = x_eq + (efrac - x_eq) * (1 - exp(-t_frac)) / t_frac;

	  /* Correct electron and hydrogen species for this estimate */

	  total_h = BaryonField[HINum][index] + BaryonField[HIINum][index];
	  new_hii = x_eq * total_h;

//	  printf("(%"ISYM" %"ISYM" %"ISYM") t_i=%.2g, x_eq=%.2g, t_frac=%.2g, "
//		 "x_est=%.2g, x0=%.2g, kph=%.2g\n",
//		 i, j, k, t_i, x_eq, t_frac, x_estimate, efrac, 
//		 BaryonField[kphHINum][index]);

	  BaryonField[DeNum][index] += new_hii - BaryonField[HIINum][index];
	  BaryonField[HINum][index] = (1.0 - x_eq) * total_h;
	  BaryonField[HIINum][index] = new_hii;

	} // ENDIF shell cell
	
      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}
