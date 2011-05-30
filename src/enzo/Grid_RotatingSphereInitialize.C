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


  /* declarations */
 
  int size = 1, dim, cellindex;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  FLOAT r,x,y,z, radius, zdist;

  float sintheta, costheta, omega;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

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
	  POW(y-RotatingSphereCenterPosition[1], 2.0);

	radius = sqrt(radius);  // ok, now it's just radius

	zdist = fabs(z-RotatingSphereCenterPosition[2]);

	if ( (radius <= RotatingSphereCoreRadius) && (zdist <= RotatingSphereCoreRadius) ){

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
	  
	} // if (r <= RotatingSphereCoreRadius)

      } // for (i = 0; i < GridDimension[0]; i++)

  if(debug){

    printf("Exiting RotatingSphereInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}

