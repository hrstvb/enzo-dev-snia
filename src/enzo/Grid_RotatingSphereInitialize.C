/***********************************************************************
/
/  GRID CLASS (INITIALIZE ROTATING SPHERE TEST)
/
/  written by: Brian O'Shea
/  date:       May 2011
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
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

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::RotatingSphereInitializeGrid(FLOAT RotatingSphereCoreRadius,
				       FLOAT RotatingSphereCenterPosition[MAX_DIMENSION],
				       float RotatingSphereLambda,
				       float RotatingSphereCentralDensity,
				       float RotatingSphereCentralTemperature)
{
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if(debug){
    printf("Entering RotatingSphereInitializeGrid\n");
    fflush(stdout);
  }
 
  printf("RotatingSphereCoreRadius = %e\n",RotatingSphereCoreRadius);
  printf("RotatingSphereCenterPosition = %e %e %e\n", 
	 RotatingSphereCenterPosition[0],
	 RotatingSphereCenterPosition[1],
	 RotatingSphereCenterPosition[2]);
  printf("RotatingSphereLambda = %e\n",RotatingSphereLambda);
  printf("RotatingSphereCentralDensity = %e\n",RotatingSphereCentralDensity);
  printf("RotatingSphereCentralTemperature = %e\n",RotatingSphereCentralTemperature);


  /* declarations */
 
  int size = 1, dim, cellindex;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  FLOAT r,x,y,z, radius, zdist;

  float sintheta, costheta, omega;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;


  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, TimeUnits=1.0,
    VelocityUnits=1.0;
  double MassUnits=1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }



  /* set fields in the cylinder region */
 
  int index, jndex, i, j, k;
  float outside_rho, outside_TE, outside_GE;

  outside_rho =  BaryonField[DensNum][0];

  // updated to include correct gravitational constant and more accurate constant (corrections by J-H Choi, U. Kentucky)
  omega = RotatingSphereLambda * sqrt((GravitationalConstant / (4.0*M_PI)) * RotatingSphereCentralDensity * outside_rho) / 0.146;

  if(HydroMethod==2){  // ZEUS

    outside_TE = BaryonField[TENum][0];

  } else { // PPM

    outside_TE = BaryonField[TENum][0];
    
    if(DualEnergyFormalism){
      outside_GE = BaryonField[GENum][0];
    }

  }  // if(HydroMethod==2)

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++){
 
	/* Compute position */
	x=y=z=0.0;

	cellindex = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	/* Find distance from center. */

	// it's REALLY r^2 right now
	radius = POW(x-RotatingSphereCenterPosition[0], 2.0) +
	  POW(y-RotatingSphereCenterPosition[1], 2.0) +
	  POW(z-RotatingSphereCenterPosition[2], 2.0);

	radius = sqrt(radius);  // ok, now it's just radius


	BaryonField[DensNum][cellindex] = outside_rho * RotatingSphereCentralDensity;

	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
	  BaryonField[MetalNum][cellindex] = BaryonField[DensNum][cellindex]*TestProblemData.MetallicityField_Fraction;

	sintheta = (y-RotatingSphereCenterPosition[1])/radius;
	costheta = (x-RotatingSphereCenterPosition[0])/radius;


	// x,y, and maybe z velocity.  
	BaryonField[Vel1Num][cellindex] = -1.0*sintheta*omega*radius;

	BaryonField[Vel2Num][cellindex] = costheta*omega*radius;

	BaryonField[Vel3Num][cellindex] = 0.0;

	if(HydroMethod == 2){

	  // ZEUS
	  BaryonField[TENum][cellindex] = outside_TE / RotatingSphereCentralDensity;

	} else {
	    
	  // PPM
	  BaryonField[TENum][cellindex] = outside_TE / RotatingSphereCentralDensity
	    + 0.5 * BaryonField[Vel1Num][cellindex] * BaryonField[Vel1Num][cellindex]
	    + 0.5 * BaryonField[Vel2Num][cellindex] * BaryonField[Vel2Num][cellindex]
	    + 0.5 * BaryonField[Vel3Num][cellindex] * BaryonField[Vel3Num][cellindex];
	    
	  // gas energy (PPM dual energy formalims)
	  if(DualEnergyFormalism)
	    BaryonField[GENum][cellindex] = outside_GE / RotatingSphereCentralDensity;
	    
	} // if(HydroMethod == 2)
	  

	/* set multispecies values --- EVERYWHERE, not just inside the sphere radius! */

	if(TestProblemData.MultiSpecies) {

	  BaryonField[HIINum][cellindex] = TestProblemData.HII_Fraction * 
	    TestProblemData.HydrogenFractionByMass * BaryonField[DensNum][cellindex];
	      
	  BaryonField[HeIINum][cellindex] = TestProblemData.HeII_Fraction *
	    BaryonField[DensNum][cellindex] * (1.0-TestProblemData.HydrogenFractionByMass);
	      
	  BaryonField[HeIIINum][cellindex] = TestProblemData.HeIII_Fraction *
	    BaryonField[DensNum][cellindex] * (1.0-TestProblemData.HydrogenFractionByMass);

	  BaryonField[HeINum][cellindex] = 
	    (1.0 - TestProblemData.HydrogenFractionByMass)*BaryonField[DensNum][cellindex] -
	    BaryonField[HeIINum][cellindex] - BaryonField[HeIIINum][cellindex];
	      
	  if(TestProblemData.MultiSpecies > 1){
	    BaryonField[HMNum][cellindex] = TestProblemData.HM_Fraction *
	      BaryonField[HIINum][cellindex];
		
	    BaryonField[H2INum][cellindex] = TestProblemData.H2I_Fraction *
	      BaryonField[0][cellindex] * TestProblemData.HydrogenFractionByMass;
		
	    BaryonField[H2IINum][cellindex] = TestProblemData.H2II_Fraction * 2.0 *
	      BaryonField[HIINum][cellindex];
	  }

	  // HI density is calculated by subtracting off the various ionized fractions
	  // from the total
	  BaryonField[HINum][cellindex] = TestProblemData.HydrogenFractionByMass*BaryonField[0][cellindex]
	    - BaryonField[HIINum][cellindex];
	  if (MultiSpecies > 1)
	    BaryonField[HINum][cellindex] -= (BaryonField[HMNum][cellindex] + BaryonField[H2IINum][cellindex]
					      + BaryonField[H2INum][cellindex]);

	  // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
	  // density for convenience) is calculated by summing up all of the ionized species.
	  // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
	  // calculating mass density, not number density (because the BaryonField values are 4x as
	  // heavy for helium for a single electron)
	  BaryonField[DeNum][cellindex] = BaryonField[HIINum][cellindex] +
	    0.25*BaryonField[HeIINum][cellindex] + 0.5*BaryonField[HeIIINum][cellindex];
	  if (MultiSpecies > 1)
	    BaryonField[DeNum][cellindex] += 0.5*BaryonField[H2IINum][cellindex] -
	      BaryonField[HMNum][cellindex];
	      
	  // Set deuterium species (assumed to be a negligible fraction of the total, so not
	  // counted in the conservation)
	  if(TestProblemData.MultiSpecies > 2){
	    BaryonField[DINum ][cellindex]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][cellindex];
	    BaryonField[DIINum][cellindex] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][cellindex];
	    BaryonField[HDINum][cellindex] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][cellindex];
	  }

	} // if(TestProblemData.MultiSpecies)	



      } // for (i = 0; i < GridDimension[0]; i++)

  if(debug){

    printf("Exiting RotatingSphereInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}

