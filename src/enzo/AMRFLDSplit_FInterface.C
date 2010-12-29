/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, Fortran interfaces
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(amrfldsplit_setupsystem)(
   Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, float *E0, float *E, float *kappaE, 
   float *eta, float *dt, FLOAT *a, FLOAT *a0, FLOAT *adot, FLOAT *adot0, 
   int *ESpectrum, float *theta, float *aUnits, float *LenUnits, 
   float *LenUnits0, float *ErUnits, float *ErUnits0, float *NiUnits, 
   float *NiUnits0, int *rank, float *dx, float *dy, float *dz, int *BCTypeXl, 
   int *BCTypeXr, int *BCTypeYl, int *BCTypeYr, int *BCTypeZl, int *BCTypeZr, 
   int *x0s, int *x0e, int *x1s, int *x1e, int *x2s, int *x2e, int *Nx, int *Ny, 
   int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, 
   int *xlface, int *xrface, int *ylface, int *yrface, int *zlface, int *zrface, 
   int *ier);

extern "C" void FORTRAN_NAME(amrfldsplit_opacity)(
   float *kappaE, float *time, float *rho, float *n_HI, float *n_HeI, 
   float *n_HeII, FLOAT *a, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, float *IsEsHeIInu, 
   float *aUnits, float *DenUnits, float *LenUnits, float *TimeUnits, 
   float *NiUnits, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(amrfldsplit_radiationsource)(
   float *Ersrc, float *time, FLOAT *a, int *ProblemType, int *ESpectrum, 
   float *NGammaDot, float *EtaRadius, float *EtaCenter, float *aUnits, 
   float *LenUnits, float *TimeUnits, float *ErUnits, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);




/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the AMRFLDSplit class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int AMRFLDSplit::SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, 
			     float *E0, float *E, float *kappaE, float *eta) 
{
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int x0s=1, x0e=LocDims[0], x1s=1, x1e=LocDims[1], x2s=1, x2e=LocDims[2];
  int ier;
  FORTRAN_NAME(amrfldsplit_setupsystem)
    (mat, rhs, rhsnorm, E0, E, kappaE, eta, &dt, &a, &a0, &adot, &adot0, 
     &ESpectrum, &theta, &aUnits, &LenUnits, &LenUnits0, &ErUnits, &ErUnits0, 
     &NiUnits, &NiUnits0, &rank, &dx[0], &dx[1], &dx[2], &BdryType[0][0], 
     &BdryType[0][1], &BdryType[1][0], &BdryType[1][1], &BdryType[2][0], 
     &BdryType[2][1], &x0s, &x0e, &x1s, &x1e, &x2s, &x2e, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &xlface, &xrface, 
     &ylface, &yrface, &zlface, &zrface, &ier);

  // combine the processor-local rhsnorm values together before returning 
  // (since I perform no MPI in F90 modules)
  float rhssum=0.0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(rhsnorm, &rhssum, 1, DataType, MPI_SUM, MPI_COMM_WORLD);
#else
  rhssum = *rhsnorm;
#endif
  *rhsnorm = sqrt(rhssum);

  return(ier);
}


/********/
int AMRFLDSplit::Opacity(float *kappaE, float *time)
{
  int ier;
  FLOAT aval = (a+a0)*0.5;
  float dUn  = (DenUnits + DenUnits0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  float nUn  = (NiUnits + NiUnits0)*0.5;
  FORTRAN_NAME(amrfldsplit_opacity)
    (kappaE, time, rho, HI, HeI, HeII, &aval, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, &intSigESigHeII, 
     &intSigESigHeIInu, &aUnits, &dUn, &lUn, &TimeUnits, &nUn, &Nchem, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int AMRFLDSplit::RadiationSource(float *Ersrc, float *time)
{
  int ier;
  FLOAT aval = (a+a0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  float rUn  = (ErUnits + ErUnits0)*0.5;
  FORTRAN_NAME(amrfldsplit_radiationsource)
    (Ersrc, time, &aval, &ProblemType, &ESpectrum, &NGammaDot, &EtaRadius, 
     EtaCenter, &aUnits, &lUn, &TimeUnits, &rUn, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], 
     &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


#endif
