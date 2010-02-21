#ifdef USE_COSMOTABLE

/***********************************************************************
/
/ CALCULATES A LOOKUP TABLE FOR 'Time*H0' and '(da/dt)/H0' AS A FUNCTION
/ OF SCALE FACTOR.
/
/  written by: Michael Kuhlen
/  date:       February 2010
/
/
/  PURPOSE:
/
/  NOTE: some routines adopted and adapted from PKDGRAV's cosmo.c.
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

// needed for dRombergO
#define EPSCOSMO 1e-5
#define MAXLEV 13
#define FLOAT_MAXVAL 1e20
#define MAX_ITER 100

//double arccosh(double x);
//double arcsinh(double x);

double csmExp2Time(double ScaleFactor);
double csmExp2Hub(double ScaleFactor);
double csmTime2Exp(double Time);
double csmTime2Hub(double Time);
double csmCosmoTint(double Y);
double dRombergO(double (*func)(double),double a,double b,double eps);


int CosmologyCalculateTable(void)
{  
  int i;

  FLOAT ScaleFactor, MinScaleFactor, MaxScaleFactor;

//   printf("ScaleFactor = %e\n",ScaleFactor);
//   printf("HubbleConstantNow = %e\n",HubbleConstantNow);
//   printf("OmegaMatterNow = %e\n",OmegaMatterNow);
//   printf("OmegaRadiationNow = %e\n",OmegaRadiationNow);
//   printf("OmegaLambdaNow = %e\n",OmegaLambdaNow);

        
  // Note: don't know where to free this memory...
  CosmologyTable_a = new FLOAT[CosmologyTableSize];
  CosmologyTable_dadt = new FLOAT[CosmologyTableSize];
  CosmologyTable_time = new FLOAT[CosmologyTableSize];
  
  MinScaleFactor = 1.0/1001.0;
  MaxScaleFactor = 1.0;
  
  CosmologyTableDelta = (MaxScaleFactor - MinScaleFactor) / (FLOAT)(CosmologyTableSize - 1);
  
  for(i=0, ScaleFactor=MinScaleFactor; i<CosmologyTableSize; i++, ScaleFactor+=CosmologyTableDelta) {
    CosmologyTable_a[i] = ScaleFactor;

    CosmologyTable_dadt[i] = csmExp2Hub(ScaleFactor) * ScaleFactor / HubbleConstantNow;      
    
    CosmologyTable_time[i] = csmExp2Time(ScaleFactor) * HubbleConstantNow;    
  }
  
#ifdef USE_MPI
  printf("MyProcessorNumber = %d\n",MyProcessorNumber);
  if (CosmologyTableWriteToFile && MyProcessorNumber == ROOT_PROCESSOR) {
#else
  if (CosmologyTableWriteToFile) {
#endif

    FILE *fptr;
    if ((fptr = fopen(CosmologyTableFilename, "w")) == NULL) {
      fprintf(stderr, "Error opening CosmologyTable output file %s\n", CosmologyTableFilename);
      ENZO_FAIL("");
    }
    
    fprintf(fptr,"#          a           z    (da/dt)/H0       H0*Time\n");
    fprintf(fptr,"#\n");
    
    for(i=0;i<CosmologyTableSize;i++) {
      FLOAT Redshift = 1.0/CosmologyTable_a[i] - 1.0;
      fprintf(fptr,"%12.6f %11.6f  %12.6e  %12.6e\n",CosmologyTable_a[i],Redshift,CosmologyTable_dadt[i],CosmologyTable_time[i]);
    }
    
    fclose(fptr);
  }
  
  CosmologyTableCalculated = TRUE;

  return SUCCESS;
}


double csmTime2Exp(double Time) {
  double al=0,ah=1,a0,a1=1,at,a;
  double th,f,f1,h,ho;
  int j;
  
  assert(Time > 0);
  th = csmExp2Time(ah);
  /*
  ** Search for upper bracket if needed.
  */
  while (Time > th) {
    a0 = a1;
    a1 = ah;
    ah = a1+a0;
    th = csmExp2Time(ah);
  }
  a = 0.5*(al+ah);
  ho = ah-al;
  h = ho;
  f = Time - dRombergO(csmCosmoTint,0.0,POW(a,1.5),EPSCOSMO);
  f1 = 1/(a*csmExp2Hub(a));
  for (j=0;j<MAX_ITER;++j) {
    if (a+f/f1 < al || a+f/f1 > ah || fabs(2*f) > fabs(ho*f1)) {
      /*
      ** Bisection Step.
      */
      ho = h;
      h = 0.5*(ah-al);
      a = al+h;
      /*
	printf("bisect al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
      */
      if (a == al) return a;
    }
    else {
      /*
      ** Newton Step.
      */
      ho = h;
      h = f/f1;
      at = a;
      a += h;
      /*
	printf("newton al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
      */
      if (a == at) return a;
    }
    if (fabs(h) < EPSCOSMO) {
      /*
	printf("converged al:%.14g ah:%.14g a:%.14g t:%.14g == %.14g\n",
	al,ah,a,dRombergO(csmCosmoTint,0.0,POW(a,1.5),EPSCOSMO*1e-1),
	Time);
      */
      return a;
    }
    f = Time - dRombergO(csmCosmoTint,0.0,POW(a,1.5),EPSCOSMO*1e-1);
    f1 = 1/(a*csmExp2Hub(a));
    if (f < 0) ah = a;
    else al = a;
  }
  assert(0);
  
  return 0.0; /* We never reach here, but this keeps the compiler happy */
}

double csmExp2Time(double ScaleFactor) {
  double Time;

  Time = dRombergO(csmCosmoTint, 0.0, POW(ScaleFactor, 1.5), EPSCOSMO);
  return Time;
}

double csmCosmoTint(double Y) {
  double ScaleFactor = POW(Y, 2.0/3.0);
  
  assert(ScaleFactor > 0.0);
  return 2.0/(3.0*Y*csmExp2Hub(ScaleFactor));
}

double csmExp2Hub(double ScaleFactor) {
  double OmegaCurvatureNow = 1.0 - OmegaMatterNow - OmegaLambdaNow - OmegaRadiationNow;
  double a2,a3,a4;
  double val;

  assert(ScaleFactor > 0.0);
  
  a2 = ScaleFactor*ScaleFactor;
  a3 = a2 * ScaleFactor;
  a4 = a3 * ScaleFactor;

  //  printf("ScaleFactor = %e\n",ScaleFactor);

  val = HubbleConstantNow * sqrt( OmegaMatterNow/a3 +
				  OmegaCurvatureNow/a2 +
				  OmegaRadiationNow/a4 + 
				  OmegaLambdaNow );
  //  printf("val = %e\n",val);

  return val;

}


double csmTime2Hub(double Time) {
  double ScaleFactor = csmTime2Exp(Time) * (1 + InitialRedshift);
  
  assert(ScaleFactor > 0.0);
  return csmExp2Hub(ScaleFactor);
}


/*
** Romberg integrator for an open interval.
*/

double dRombergO(double (*func)(double),double a,double b,double eps) {
  double tllnew;
  double tll;
  double tlk[MAXLEV+1];
  Eint32 n = 1;
  Eint32 nsamples = 1;
  
  tlk[0] = tllnew = (b-a)*(*func)(0.5*(b+a));
  if (a == b) return tllnew;
  tll = FLOAT_MAXVAL;
  
  eps*=0.5;
  
  while ((fabs((tllnew-tll)/(tllnew+tll)) > eps) && (n < MAXLEV)) {
    /*
     * midpoint rule.
     */
    double deltax;
    double tlktmp;
    Eint32 i;
    
    nsamples *= 3;
    deltax = (b-a)/nsamples;
    tlktmp = tlk[0];
    tlk[0] = tlk[0]/3.0;
    
    for (i=0;i<nsamples/3;i++) {
      tlk[0] += deltax*(*func)(a + (3*i + 0.5)*deltax);
      tlk[0] += deltax*(*func)(a + (3*i + 2.5)*deltax);
    }
    
    /*
     * Romberg extrapolation.
     */
    
    for (i=0;i<n;i++) {
      double tlknew = (POW(9.0, i+1.)*tlk[i] - tlktmp)
	/(POW(9.0, i+1.) - 1.0);
      
      tlktmp = tlk[i+1];
      tlk[i+1] = tlknew;
    }
    tll = tllnew;
    tllnew = tlk[n];
    n++;
  }
  
  assert((fabs((tllnew-tll)/(tllnew+tll)) < eps));
  
  return tllnew;
}


#endif
