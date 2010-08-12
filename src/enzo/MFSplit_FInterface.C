/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Frequency, Multi-species, Split Problem Class
/  Fortran interfaces.
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef TRANSFER
#include "MFSplit.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(mfsplit_setupsystem)(
   Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, int *freq, float *E, float *Eold, 
   float *HI, float *HIold, float *HeI, float *HeIold, float *HeII, 
   float *HeIIold, float *src, int *LimType, float *dt, FLOAT *a, FLOAT *a0, 
   float *theta, float *aUn, float *lUn, float *lUn0, float *nUn, float *nUn0, 
   int *rank, float *dx, float *dy, float *dz, int *BCxl, int *BCxr, int *BCyl, 
   int *BCyr, int *BCzl, int *BCzr, int *x0s, int *x0e, int *x1s, int *x1e, 
   int *x2s, int *x2e, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *xlface, 
   int *xrface, int *ylface, int *yrface, int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_radiationsource)(
   float *eta1, float *eta2, float *eta3, float *FS_NGammaDot, float *time, 
   FLOAT *a, int *Model, int *ProblemType, int *ESpectrum, int *Nchem, 
   int *SMEmiss, float *NGammaDot, float *chibar, float *chinuint, 
   float *EtaRadius, float *EtaCenter, float *aUn, float *dUn, float *vUn, 
   float *lUn, float *tUn, float *E1Un, float *E2Un, float *E3Un, float *eUn, 
   float *nUn, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, 
   float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_gasenergysource)(
   float *ecsrc, float *time, FLOAT *a, int *Model, int *ProblemType, 
   int *Nchem, float *aUn, float *dUn, float *vUn, float *lUn, 
   float *tUn, float *E1Un, float *E2Un, float *E3Un, float *eUn, 
   float *nUn, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, 
   float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_chemistrysource)(
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *time, FLOAT *a, 
   int *Model, int *ProblemType, int *Nchem, float *aUn, float *dUn, 
   float *vUn, float *lUn, float *tUn, float *E1Un, float *E2Un, float *E3Un, 
   float *eUn, float *nUn, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, 
   float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_radinit)(
   float *Ef, float *E1, float *E2, float *E3, float *E1Units, 
   float *E2Units, float *E3Units, float *chiint, float *chinuint, 
   int *ESpectrum, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_computeradiationintegrals)(
   float *piHI, float *piHeI, float *piHeII, float *GHI, float *GHeI, 
   float *GHeII, float *Ef, float *E1, float *E2, float *E3, float *E1old, 
   float *E2old, float *E3old, int *Nchem, int *ESpectrum, float *chibar, 
   float *fsUn, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_enforceradiationbounds)(
   float *Ef, float *E1, float *E2, float *E3, float *fsUn, float *E1Un, 
   float *E2Un, float *E3Un, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_analyticchemistry)(
   float *ec, float *HI, float *HeI, float *HeII, float *HIold, 
   float *HeIold, float *HeIIold, float *dt, float *vx, float *vy, 
   float *vz, float *rho, float *eh, float *ecsrc, float *HIsrc, 
   float *HeIsrc, float *HeIIsrc, float *piHI, float *piHeI, 
   float *piHeII, float *GHI, float *GHeI, float *GHeII, float *gamma, 
   float *HFrac, int *model, int *DualEnergy, FLOAT *a, FLOAT *a0, 
   FLOAT *adot, FLOAT *adot0, float *CompA, float *Comp_xray, 
   float *Comp_temp, int *NTempBins, float *TStart, float *TEnd, 
   float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, 
   float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, 
   float *ciHITb, float *ciHeITb, float *ciHeISTb, float *ciHeIITb, 
   float *reHIITb, float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, 
   float *bremTb, float *aUn, float *dUn, float *dUn0, float *vUn, 
   float *lUn, float *lUn0, float *eUn, float *nUn, float *nUn0, 
   float *ecScale, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_analyticinitguess)(
   float *Ef, float *E1, float *E2, float *E3, float *ec, float *HI, 
   float *HeI, float *HeII, float *dt, float *vx, float *vy, float *vz, 
   float *rho, float *eh, float *E1src, float *E2src, float *E3src, 
   float *ecsrc, float *HIsrc, float *HeIsrc, float *HeIIsrc, float *gamma, 
   float *hfrac, int *model, int *dualenergy, int *ESpectrum, FLOAT *a, 
   FLOAT *a0, FLOAT *adot, FLOAT *adot0, float *CompA, float *Comp_xray, 
   float *Comp_temp, float *piHI, float *piHeI, float *piHeII, float *GHI, 
   float *GHeI, float *GHeII, int *NTempBins, float *TStart, float *TEnd, 
   float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, 
   float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, 
   float *ciHeITb, float *ciHeISTb, float *ciHeIITb, float *reHIITb, 
   float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, float *bremTb, 
   float *chibar, float *aUn, float *dUn, float *dUn0, float *vUn, float *lUn, 
   float *lUn0, float *fsUn, float *eUn, float *nUn, float *nUn0, 
   float *ecScale, int *Nchem, float *dx, float *dy, float *dz, int *Nx, 
   int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, 
   int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfsplit_fsradiationmask)(
   float *Ef, float *E1, float *E2, float *E3, float *fsUn, float *E1Un, 
   float *E2Un, float *E3Un, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);




/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the MFSplit class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int MFSplit::SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, 
			 int freq, EnzoVector *u, EnzoVector *srcvec,
			 float thisdt) 
{
  float stime = MPI_Wtime();
  int x0s=1, x0e=LocDims[0], x1s=1, x1e=LocDims[1], x2s=1, x2e=LocDims[2];
  x0s -= (OnBdry[0][0] && (BdryType[0][0]==1)) ? 1 : 0;
  x0e += (OnBdry[0][1] && (BdryType[0][1]==1)) ? 1 : 0;
  x1s -= (OnBdry[1][0] && (BdryType[1][0]==1)) ? 1 : 0;
  x1e += (OnBdry[1][1] && (BdryType[1][1]==1)) ? 1 : 0;
  x2s -= (OnBdry[2][0] && (BdryType[2][0]==1)) ? 1 : 0;
  x2e += (OnBdry[2][1] && (BdryType[2][1]==1)) ? 1 : 0;
  int ier;
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  float *E, *Eold, *src;
  if (freq == 1) {
    E = u->GetData(iE1);
    Eold = U0->GetData(iE1);
    src = srcvec->GetData(iE1);
  }
  else if (freq == 2) {
    E = u->GetData(iE2);
    Eold = U0->GetData(iE2);
    src = srcvec->GetData(iE2);
  }
  else if (freq == 3) {
    E = u->GetData(iE3);
    Eold = U0->GetData(iE3);
    src = srcvec->GetData(iE3);
  }
  else
    ENZO_VFAIL("MFSplit FInterface: freq = %"ISYM" undefined!\n",freq)
  FORTRAN_NAME(mfsplit_setupsystem)
    (mat, rhs, rhsnorm, &freq, E, Eold, HI, HIold, HeI, HeIold, HeII, HeIIold, 
     src, &LimType, &thisdt, &a, &a0, &theta, &aUnits, &LenUnits, &LenUnits0, 
     &nUnits, &nUnits0, &rank, &dx[0], &dx[1], &dx[2], &(BdryType[0][0]), 
     &(BdryType[0][1]), &(BdryType[1][0]), &(BdryType[1][1]), &(BdryType[2][0]), 
     &(BdryType[2][1]), &x0s, &x0e, &x1s, &x1e, &x2s, &x2e, &Nchem, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &xlface, &xrface, &ylface, 
     &yrface, &zlface, &zrface, &ier);
  float ftime = MPI_Wtime();
  timers[0] += ftime - stime;
  return(ier);
}

/********/
int MFSplit::RadiationSource(EnzoVector *extsrc, float *time, 
			     float *FS_NGammaDot)
{
  float stime = MPI_Wtime();
  int ier;
  float *eta1 = extsrc->GetData(iE1);
  float *eta2 = extsrc->GetData(iE2);
  float *eta3 = extsrc->GetData(iE3);
  FORTRAN_NAME(mfsplit_radiationsource)
    (eta1, eta2, eta3, FS_NGammaDot, time, &a, &Model, &ProblemType, 
     &ESpectrum, &Nchem, &StarMakerEmissivityField, &NGammaDot, &chibar, 
     &chinuint, &EtaRadius, &EtaCenter[0], &aUnits, &DenUnits, &VelUnits, 
     &LenUnits, &TimeUnits, &E1Units, &E2Units, &E3Units, &eUnits, &nUnits, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[5] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::GasEnergySource(EnzoVector *extsrc, float *time)
{
  float stime = MPI_Wtime();
  int ier;
  float *ecsrc = extsrc->GetData(iec);
  FORTRAN_NAME(mfsplit_gasenergysource)
    (ecsrc, time, &a, &Model, &ProblemType, &Nchem, &aUnits, &DenUnits, 
     &VelUnits, &LenUnits, &TimeUnits, &E1Units, &E2Units, &E3Units, &eUnits, 
     &nUnits, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], 
     &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[5] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::ChemistrySource(EnzoVector *extsrc, float *time)
{
  float stime = MPI_Wtime();
  int ier;
  float *HIsrc = extsrc->GetData(iHI);
  float *HeIsrc = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  FORTRAN_NAME(mfsplit_chemistrysource)
    (HIsrc, HeIsrc, HeIIsrc, time, &a, &Model, &ProblemType, &Nchem, &aUnits, 
     &DenUnits, &VelUnits, &LenUnits, &TimeUnits, &E1Units, &E2Units, &E3Units, 
     &eUnits, &nUnits, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[5] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::RadInit(float *Ef, float *E1, float *E2, float *E3, 
		     float *E1Units, float *E2Units, float *E3Units, 
		     float *chiint, float *chinuint)
{
  int ier;
  FORTRAN_NAME(mfsplit_radinit)(Ef, E1, E2, E3, E1Units, E2Units, E3Units, 
				chiint, chinuint, &ESpectrum, &ier);
  return(ier);
}


/********/
int MFSplit::ComputeRadiationIntegrals(EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *E1old = U0->GetData(iE1);
  float *E2old = U0->GetData(iE2);
  float *E3old = U0->GetData(iE3);

  FORTRAN_NAME(mfsplit_computeradiationintegrals)
    (piHI, piHeI, piHeII, GHI, GHeI, GHeII, Efree, E1, E2, E3, E1old, E2old, 
     E3old, &Nchem, &ESpectrum, &chibar, &fsUnits, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[6] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::EnforceRadiationBounds(float *Ef, float *E1, float *E2, float *E3)
{
  float stime = MPI_Wtime();
  int ier;
  FORTRAN_NAME(mfsplit_enforceradiationbounds)
    (Ef, E1, E2, E3, &fsUnits, &E1Units, &E2Units, &E3Units, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[7] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::AnalyticChemistry(EnzoVector *sol, EnzoVector *src, float thisdt)
{
  float stime = MPI_Wtime();
  int ier;
  float *ec = sol->GetData(iec);
  float *HI = sol->GetData(iHI);
  float *HeI = sol->GetData(iHeI);
  float *HeII = sol->GetData(iHeII);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  float *ecsrc = src->GetData(iec);
  float *HIsrc = src->GetData(iHI);
  float *HeIsrc = src->GetData(iHeI);
  float *HeIIsrc = src->GetData(iHeII);
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  FORTRAN_NAME(mfsplit_analyticchemistry)
    (ec, HI, HeI, HeII, HIold, HeIold, HeIIold, &thisdt, vx, vy, vz, rho, 
     eh, ecsrc, HIsrc, HeIsrc, HeIIsrc, piHI, piHeI, piHeII, GHI, GHeI, GHeII, 
     &Gamma, &HFrac, &Model, &dualenergy, &a, &a0, &adot, &adot0, 
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &aUnits, &DenUnits, &DenUnits0, 
     &VelUnits, &LenUnits, &LenUnits0, &eUnits, &nUnits, &nUnits0, &eScale, 
     &Nchem, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &ier);
  float ftime = MPI_Wtime();
  timers[8] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::AnalyticInitGuess(EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *ec = u->GetData(iec);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *E1src = extsrc->GetData(iE1);
  float *E2src = extsrc->GetData(iE2);
  float *E3src = extsrc->GetData(iE3);
  float *ecsrc = extsrc->GetData(iec);
  float *HIsrc = extsrc->GetData(iHI);
  float *HeIsrc = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  FORTRAN_NAME(mfsplit_analyticinitguess)
    (Efree, E1, E2, E3, ec, HI, HeI, HeII, &dt, vx, vy, vz, rho, eh, 
     E1src, E2src, E3src, ecsrc, HIsrc, HeIsrc, HeIIsrc, &Gamma, 
     &HFrac, &Model, &dualenergy, &ESpectrum, &a, &a0, &adot, &adot0, 
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, piHI, piHeI, 
     piHeII, GHI, GHeI, GHeII, &CoolData.NumberOfTemperatureBins, 
     &CoolData.TemperatureStart, &CoolData.TemperatureEnd, RateData.k1, 
     RateData.k2, RateData.k3, RateData.k4, RateData.k5, RateData.k6, 
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI, 
     CoolData.ciHeI, CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
     CoolData.reHeII1, CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, 
     &chibar, &aUnits, &DenUnits, &DenUnits0, &VelUnits, &LenUnits, &LenUnits0, 
     &fsUnits, &eUnits, &nUnits, &nUnits0, &eScale, &Nchem, &dx[0], &dx[1], 
     &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[9] += ftime - stime;
  return(ier);
}


/********/
int MFSplit::FSRadiationMask(float *Ef, float *E1, float *E2, float *E3)
{
  int ier;
  FORTRAN_NAME(mfsplit_fsradiationmask)
    (Ef, E1, E2, E3, &fsUnits, &E1Units, &E2Units, &E3Units, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


#endif
