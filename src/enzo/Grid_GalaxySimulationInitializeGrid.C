/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/  modified1:  Elizabeth Tasker, Oct, 2006 (tidied up)
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
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

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]

int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, double *MassUnits, FLOAT Time);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* Internal routines */

float gauss_mass(FLOAT rCell, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float CircularVelocity, float SoundSpeed, FLOAT rCore, FLOAT ScaleHeightz0, FLOAT cellwidth);
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);

static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits;
static double MassUnits;

int grid::GalaxySimulationInitializeGrid(FLOAT DiskRadius,
					 FLOAT ExternalGravityRadius,
					 float ExternalGravityOrientation[MAX_DIMENSION],
					 float ExternalGravityConstant,
					 FLOAT DiskPosition[MAX_DIMENSION], 
					 FLOAT ScaleHeightz,
					 float DiskTemperature,
					 float InitialTemperature,
					 float AngularMomentum[MAX_DIMENSION],
					 float UniformVelocity[MAX_DIMENSION], 
					 int UseMetallicityField, 
					 float GalaxySimulationInflowTime,
					 float GalaxySimulationInflowDensity,
					 int level)
{
 /* declarations */

 int dim, i, j, k, m, field, disk, size, MetalNum, vel;
 int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
   DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
 float DiskDensity, CircularVelocity, DiskVelocityMag, SoundSpeed;
 FLOAT rCore;
  
  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (HydroMethod == MHD_RK) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
    FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    if (UseDivergenceCleaning) {
      FieldType[NumberOfBaryonFields++] = Phi_pField;
    }
  }

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  if (UseMetallicityField)
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* fake it with metals */

 /* Return if this doesn't concern us. */

 if (ProcessorNumber != MyProcessorNumber) 
   return SUCCESS;


 /* Set various units. */

 float CriticalDensity = 1, BoxLength = 1, mu = 0.6;
 FLOAT a, dadt, ExpansionFactor = 1;
 if (ComovingCoordinates) {
   CosmologyComputeExpansionFactor(Time, &a, &dadt);
   ExpansionFactor = a/(1.0+InitialRedshift);
   CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);
   CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
   BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
 } else {
   GetUnits(&DensityUnits,&LengthUnits,
            &TemperatureUnits,&TimeUnits,
            &VelocityUnits, &MassUnits, Time);
   BoxLength = LengthUnits/Mpc;
 }

 /* Set up inflow */
 if (GalaxySimulationInflowTime > 0.0){
   TimeActionType[0] = 2;
   TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
   TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
 }

 /* compute size of fields */

 size = 1;
 for (dim = 0; dim < GridRank; dim++)
   size *= GridDimension[dim];

 /* allocate fields */

 for (field = 0; field < NumberOfBaryonFields; field++)
   if (BaryonField[field] == NULL)
     BaryonField[field] = new float[size];

 /* Loop over the mesh. */

 float density, dens1, Velocity[MAX_DIMENSION],
   temperature, temp1;
 FLOAT r, x, y = 0, z = 0;
 int n = 0;

 for (k = 0; k < GridDimension[2]; k++)
   for (j = 0; j < GridDimension[1]; j++)
     for (i = 0; i < GridDimension[0]; i++, n++) {

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	if (ComovingCoordinates)
	  density = 1.0; 
	else 
	  density = 1e-4; // A small density (in density units)
	
	temperature = temp1 = InitialTemperature;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	/* Find distance from center. */

	r = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		 POW(fabs(y-DiskPosition[1]), 2) +
		 POW(fabs(z-DiskPosition[2]), 2) );
	r = max(r, 0.1*CellWidth[0][0]);

	if (r < DiskRadius) {

	  FLOAT xpos, ypos, zpos, zheight, drad; 
	  float CellMass;
	  FLOAT xhat[3];
	  FLOAT yhat[3];

	  /* Loop over dims if using Zeus (since vel's face-centered). */

	  for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
	       dim++) {

	    /* Compute position. */

	    xpos = x-DiskPosition[0] - 
	      (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	    ypos = y-DiskPosition[1] -
	      (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	    zpos = z-DiskPosition[2] -
	      (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);
	    
	    /* Compute z and r_perp (AngularMomentum is angular momentum 
	       and must have unit length). */    

	    /* magnitude of z = r.L in L direction */

	    zheight = AngularMomentum[0]*xpos + 
	              AngularMomentum[1]*ypos +
	              AngularMomentum[2]*zpos;

	    /* position in plane of disk */

	    xhat[0] = xpos - zheight*AngularMomentum[0];
	    xhat[1] = ypos - zheight*AngularMomentum[1];
	    xhat[2] = zpos - zheight*AngularMomentum[2];
	    drad = sqrt(xhat[0]*xhat[0] + xhat[1]*xhat[1] + xhat[2]*xhat[2]);


	    /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	    xhat[0] = xhat[0]/drad;
	    xhat[1] = xhat[1]/drad;
	    xhat[2] = xhat[2]/drad;

	    /* Find another vector perpendicular to r_perp and AngularMomentum */

	    yhat[0] = AngularMomentum[1]*xhat[2] - AngularMomentum[2]*xhat[1];
	    yhat[1] = AngularMomentum[2]*xhat[0] - AngularMomentum[0]*xhat[2];
	    yhat[2] = AngularMomentum[0]*xhat[1] - AngularMomentum[1]*xhat[0];

	    /* generate rotation matrix */
	    FLOAT inv[3][3],temp;
	    int i,j;
	    
	    // matrix of basis vectors in coordinate system defined by the galaxy
	    inv[0][0] = xhat[0]; inv[0][1] = yhat[0]; inv[0][2] = AngularMomentum[0];
	    inv[1][0] = xhat[1]; inv[1][1] = yhat[1]; inv[1][2] = AngularMomentum[1];
	    inv[2][0] = xhat[2]; inv[2][1] = yhat[2]; inv[2][2] = AngularMomentum[2];
	    
	    // Matrix is orthogonal by construction so inverse = transpose
	    for (i=0;i<3;i++)
	      for (j=i+1;j<3;j++)
		{
		  temp = inv[i][j];
		  inv[i][j] = inv[j][i];
		  inv[j][i] = temp;
		}

	    /* If we're above the disk, then exit. */
	    // DiskDensity = (GasMass*SolarMass/(8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))/DensityUnits;   //Code units (rho_0)

	    // 200 km s^-1 in TT09
	    CircularVelocity = ExternalGravityConstant*VelocityUnits; // [cgs]
	    // 0.5 kpc in TT09
	    rCore = ExternalGravityRadius*LengthUnits; // [cgs]
	    // 6 km s^-1 in TT09  
	    SoundSpeed = 6e5; // [cgs]

	    DiskVelocityMag = CircularVelocity*drad*LengthUnits/sqrt(pow(rCore,2) + pow((drad*LengthUnits),2))/VelocityUnits;

	    if (dim == 0)
	      {
		CellMass = gauss_mass(drad*LengthUnits,zheight*LengthUnits, xpos*LengthUnits, ypos*LengthUnits, zpos*LengthUnits, inv, 
				      CircularVelocity, SoundSpeed, rCore, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits);
		
		dens1 = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;
	      }

	    if (dens1 < density)
	      break;

	    /* Compute velocity magnitude (divided by drad). 
	       This assumes PointSourceGravityPosition and Disk center 
	       are the same. */

	    /* Compute velocty: L x r_perp. */

	    if (dim == 0 || dim == 1)
	      Velocity[0] = DiskVelocityMag*(AngularMomentum[1]*xhat[2] -
					     AngularMomentum[2]*xhat[1]);
	    if (dim == 0 || dim == 2)
	      Velocity[1] = DiskVelocityMag*(AngularMomentum[2]*xhat[0] -
					     AngularMomentum[0]*xhat[2]);
	    if (dim == 0 || dim == 3)
	      Velocity[2] = DiskVelocityMag*(AngularMomentum[0]*xhat[1] -
					     AngularMomentum[1]*xhat[0]);
	    
	  } // end: loop over dims

	   	    
	    /* If the density is larger than the background (or the previous
	       disk), then set the velocity. */
	  if (dens1 > density) {
	    density = dens1;
	    if (temp1 == InitialTemperature)
	      temp1 = DiskTemperature;
	    temperature = temp1;
	  }

	} // end: if (r < DiskRadius)
	
	/* Set density. */

	BaryonField[0][n] = density;
	
	if (UseMetallicityField)
	  for (i = 0; i < size; i++)
	    BaryonField[MetalNum][i] = 1.0e-10;

	/* Set Velocities. */

	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits/
                           ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];
	
	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);

	if (BaryonField[1][n] <= 0)
	  printf("n = %d  temp = %g   e = %g\n", 0, temperature, 
	       BaryonField[1][0]);


     } // end loop over grid

 return SUCCESS;

}

// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature
float gauss_mass(FLOAT rCell, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float CircularVelocity, float SoundSpeed, FLOAT rCore, FLOAT ScaleHeightz0, FLOAT cellwidth)
{
  
  FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
  float Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
  float xResult [5],yResult [5];
  float Mass, ToomreQ, rho0, kappa;
  FLOAT xrot,yrot,zrot,rd, ScaleHeightz;
  FLOAT xEval,yEval,zEval;
  int i,j,k;

  const FLOAT Rsun = 0.0085*Mpc;
  // Converting from full width half magnitude to scale height
  const FLOAT h0 = ScaleHeightz0;
  const FLOAT R0 = 0.0095*Mpc;

  Mass = ToomreQ = rho0 = kappa = 0;
  xrot = yrot = zrot = xEval = yEval = zEval = rd = 0;

  if (rCell < 0.002*Mpc || rCell > 0.010*Mpc)
    ToomreQ = 20.0;
  else
    ToomreQ = 1.0;

  

  for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      xEval = xpos + EvaluationPoints[i]*cellwidth/2.0;
      for (j=0;j<5;j++)
	{
	  yResult[j] = 0.0;
	  yEval = ypos + EvaluationPoints[j]*cellwidth/2.0;
	  for (k=0;k<5;k++)
	    {
	      zEval = zpos + EvaluationPoints[k]*cellwidth/2.0;

	      xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
	      yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
	      zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];

	      rd = sqrt(POW(xrot,2.0)+POW(yrot,2.0));

	      // Fitting formula for MW disk flaring from Kalberla & Kerp (2009), ARAA, Section 3.1.4
	      //ScaleHeightz = h0*exp((rd - Rsun)/R0);
	      ScaleHeightz = h0;

	      kappa = sqrt(2.)*CircularVelocity*sqrt(POW(float(rd),2.0)+2.0*POW(float(rCore),2.0))/(POW(float(rd),2.0)+POW(float(rCore),2.0));
	      rho0 = kappa*SoundSpeed/(2.0*pi*GravConst*ToomreQ*float(ScaleHeightz));

	      yResult[j] += float(cellwidth)/2.0*Weights[k]*rho0/POW(cosh(float(zrot/ScaleHeightz)),2);
	    }
	  xResult[i] += float(cellwidth)/2.0*Weights[j]*yResult[j];
	}
      Mass += float(cellwidth)/2.0*Weights[i]*xResult[i];
    }
  return Mass;
}
