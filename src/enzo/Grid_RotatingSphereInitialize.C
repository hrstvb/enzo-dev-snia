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

// routines needed to get hydrostatic equilibrium set up
static void get_dens_temp(void);
static double enclosed_mass(double r);
static double dTdr(double r, double T);
static double dens_at_r(double r);
static double drhodr(double r);
static double accel_at_r(double r);

// some global (to this file) variables needed for HSE
static double cent_dens, cent_temp, core_radius, exterior_density, maxrad, r_ambient, *rad, *nofr, *Tofr;
static int ncells;

#define DEFAULT_MU 1.22   // assume everything's neutral!

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
 
  if(debug){
    printf("RotatingSphereCoreRadius = %e\n",RotatingSphereCoreRadius);
    printf("RotatingSphereCenterPosition = %e %e %e\n", 
	   RotatingSphereCenterPosition[0],
	   RotatingSphereCenterPosition[1],
	   RotatingSphereCenterPosition[2]);
    printf("RotatingSphereLambda = %e\n",RotatingSphereLambda);
    printf("RotatingSphereCentralDensity = %e\n",RotatingSphereCentralDensity);
    printf("RotatingSphereCentralTemperature = %e\n",RotatingSphereCentralTemperature);
  }

  /* declarations */
 
  int size = 1, dim, cellindex, distindex;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  FLOAT r,x,y,z, radius, zdist;

  float sintheta, costheta, omega;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum;

  // get physical quantities, metallicity, species information
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


  // get units
  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, TimeUnits=1.0,
    VelocityUnits=1.0;
  double MassUnits=1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }


  /* calculate density, temperature profiles */

  // convert everything to CGS
  cent_dens = (double) RotatingSphereCentralDensity;
  cent_temp = (double) RotatingSphereCentralTemperature;
  core_radius = (double) RotatingSphereCoreRadius;
  exterior_density = (double) BaryonField[DensNum][0];

  cent_dens *= 1.22 * (double) DensityUnits;  // now in CGS
  exterior_density *= 1.22 * (double) DensityUnits; // now in CGS
  core_radius *= (double) LengthUnits;  // now in CGS
  maxrad = (double) LengthUnits;

  // arrays we store n(r), T(r) in 
  ncells = 1024;
  rad = new double[ncells];
  nofr = new double[ncells];
  Tofr = new double[ncells];

  // actually do the heavy lifting here
  get_dens_temp();

  // convert arrays from CGS into Enzo internal units.
  for(int i=0; i<ncells; i++){
    rad[i] /= LengthUnits;  // convert to enzo distance
    nofr[i] *= DensityUnits;  // convert to enzo-unit density (from CGS)
    Tofr[i] /= (TemperatureUnits*(Gamma-1.0)*DEFAULT_MU);  // convert from temp to internal energy
  }

  r_ambient /= LengthUnits;  // radius of sphere that we modify in Enzo distance units

  /* set fields in the sphere region */
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

	for(int ii=0;ii<ncells-2;ii++){
	  if(radius > rad[ii] && radius <= rad[ii+1] ){
	    distindex=ii;
	  }
	}
	if(radius >= rad[ncells-1])
	  distindex=ncells-1;

	if(radius <= r_ambient){
	  BaryonField[DensNum][cellindex] = nofr[distindex];
	} 
	

	if(radius <= r_ambient){
	  sintheta = (y-RotatingSphereCenterPosition[1])/radius;
	  costheta = (x-RotatingSphereCenterPosition[0])/radius;

	  // x,y, and maybe z velocity.  
	  BaryonField[Vel1Num][cellindex] = -1.0*sintheta*omega*radius;

	  BaryonField[Vel2Num][cellindex] = costheta*omega*radius;
	  
	  BaryonField[Vel3Num][cellindex] = 0.0;
	}

	if(HydroMethod == 2){

	  BaryonField[TENum][cellindex] = Tofr[distindex];

	} else {
	    
	  // PPM
	  BaryonField[TENum][cellindex] = Tofr[distindex]
	    + 0.5 * BaryonField[Vel1Num][cellindex] * BaryonField[Vel1Num][cellindex]
	    + 0.5 * BaryonField[Vel2Num][cellindex] * BaryonField[Vel2Num][cellindex]
	    + 0.5 * BaryonField[Vel3Num][cellindex] * BaryonField[Vel3Num][cellindex];
	    
	  // gas energy (PPM dual energy formalims)
	  if(DualEnergyFormalism)
	    BaryonField[GENum][cellindex] = Tofr[distindex];
	    
	} // if(HydroMethod == 2)
	  
	// set metallicity
	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
	  BaryonField[MetalNum][cellindex] = BaryonField[DensNum][cellindex]*TestProblemData.MetallicityField_Fraction;

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

  delete [] rad;
  delete [] nofr;
  delete [] Tofr;

  return SUCCESS;

}

/* -------------------------- helper functions --------------------------- */
static void get_dens_temp(void){
  int i;
  double r_min, r_max, logdr, thisrad, thislograd, thisdens;
  r_max = maxrad;
  r_min = 1.0e-10*r_max;

  logdr = (log10(r_max) - log10(r_min))/double(ncells);

  printf("r_min, r_max, logdr, cent_dens = %e %e %e %e\n",r_min, r_max, logdr, cent_dens);

  for(int i=0; i<ncells; i++){

    thislograd = log10(r_min) + double(i)*logdr;
    thisrad = POW(10.0, thislograd);

    rad[i] = thisrad;

    if(thisrad <= core_radius){
      thisdens = cent_dens;
    } else {
      thisdens = cent_dens * POW( thisrad/core_radius, -2.0);
    }

    if(thisdens < exterior_density)
      thisdens = exterior_density;

    nofr[i] = thisdens;

    printf("i, thislograd, thisrad, thisdens, log(thisdens): %d %e %e %e %e\n", i, thislograd, rad[i], nofr[i], log10(nofr[i]));

  } //   for(int i=0; i<ncells; i++)

  fprintf(stderr,"(1)\n");

  r_ambient = -1.0;

  int amb_index = 1023;
  for(i=ncells-1; i=0; i--)
    if(nofr[i] <= exterior_density){ 
      r_ambient = rad[i];
      amb_index = i;
    }
  if(nofr[ncells-1] > exterior_density){
    amb_index=ncells-1;
    r_ambient = rad[amb_index];
  }
  fprintf(stderr,"(2) %d %e\n",amb_index, exterior_density);

  printf("r_ambient, amb_index, nofr[amb_index], rad[amb_index] = %e %d %e %e\n", r_ambient, amb_index, nofr[amb_index], rad[amb_index]);

  fprintf(stderr,"(3)\n");

  // now we do 4th order RK integration outward to calculate the temperature at any given radius
  double k1, k2, k3, k4;
  double this_T, last_T, thisdr;
  
  this_T = last_T = cent_temp;

  for(i=1; i<ncells-2; i++){

    thisdr = rad[i+1]-rad[i];

    last_T = this_T;

    // f(r,T) = dT/dr 
    
    // k1 = f(r, T)
    k1 = dTdr(rad[i], last_T);

    // k2 = f(r+0.5dr, T+0.5*k1*dr)
    k2 = dTdr(rad[i]+0.5*thisdr, last_T+0.5*k1*thisdr);

    // k3 = f(r+0.5dr, T+0.5*k2*dr)
    k3 = dTdr(rad[i]+0.5*thisdr, last_T+0.5*k2*thisdr);

    // k4 = f(r+dr, T+k3*dr)
    k4 = dTdr(rad[i]+thisdr, last_T+k3*thisdr);

    this_T = last_T + (1./6.)*thisdr*(k1 + 2.*k2 + 2.*k3 + k4);

    Tofr[i] = this_T;

    printf("i, r, n(r), T(r), Menc(r): %d %e   %e   %e   %e\n",i, rad[i],nofr[i],Tofr[i],(enclosed_mass(rad[i])/1.989e+33) );

  }

  Tofr[ncells-1] = Tofr[ncells-2];

  return;
}

/* returns enclosed mass in grams */
static double enclosed_mass(double r){
  double encmass=0.0;
  double pi = 3.14159;

  encmass = -1.0;

  if(r <= core_radius){  // inside core

    encmass = 4.0/3.0*pi*cent_dens*POW(r, 3.0); 

  } else if ( (r > core_radius) && (r <= r_ambient) ){  // where rho goes as 1/r^2

    encmass = 4.0/3.0*pi*cent_dens*POW(core_radius, 3.0)
      + 4.0*pi*cent_dens*core_radius*core_radius* (r - core_radius); 

  } else {  // outside of sphere in ambient medium

    encmass = 4.0/3.0*pi*cent_dens*POW(core_radius, 3.0)
      + 4.0*pi*cent_dens*core_radius*core_radius* (r_ambient - core_radius)
      + 4.0*pi/3.0*exterior_density*( POW(r,3.0) - POW(r_ambient,3.0) );

  }

  if(encmass<0.0){
    fprintf(stderr,"enclosed mass is < 0 %e\n",encmass);
    exit(-123);
  }

  return encmass;
}

static double dTdr(double r, double T){

  double value, bunch_of_constants;

  // 3/2 * kb / (mu*m_p)
  bunch_of_constants = 1.5 * 1.38e-16 / (1.22*1.67e-24);

  value = -dens_at_r(r) * accel_at_r(r) - bunch_of_constants*drhodr(r)*T;
  value /= (bunch_of_constants * dens_at_r(r) );

  return value;
}

static double accel_at_r(double r){

  double value;
  
  value = -6.67e-8*enclosed_mass(r)/(r*r);

  return value;

}


static double dens_at_r(double r){
  double value=-1.0; 

  if(r <= core_radius){
    value = cent_dens;
  } else if(r > core_radius && r <= r_ambient){
    value = cent_dens * POW( r/core_radius, -2.0);
  } else {
    value = exterior_density;
  }
  
  if(value < 0.0){
    fprintf(stderr,"Grid::RotatingSphereInitialize: error in dens_at_r!\n");
    exit(-123);
  }

  return value;
}

static double drhodr(double r){

  double value; 

  if(r <= core_radius){
    value = 0.0;
  } else if(r > core_radius && r <= r_ambient){
    value = -2.0*cent_dens * POW( core_radius, 2.0) * POW( r, -3.0 );
  } else {
    value = 0.0;
  }
  
  return value;
}
