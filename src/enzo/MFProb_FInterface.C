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
/  Multi-Frequency, Multi-species, Implicit Problem Class
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
#include "MFProb.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(mfprob_setupsystem)(
   Eflt64 *mat11, Eflt64 *mat12, Eflt64 *mat13, Eflt64 *mat21, Eflt64 *mat22, 
   Eflt64 *mat23, Eflt64 *mat31, Eflt64 *mat32, Eflt64 *mat33, Eflt64 *rhs1, 
   Eflt64 *rhs2, Eflt64 *rhs3, float *E1, float *E2, float *E3, float *E1old, 
   float *E2old, float *E3old, float *HI, float *HIold, float *HeI, float *HeIold, 
   float *HeII, float *HeIIold, float *adj11, float *adj12, float *adj13, 
   float *adj21, float *adj22, float *adj23, float *adj31, float *adj32, 
   float *adj33, int *LimType, int *LimImp, float *dt, float *theta, FLOAT *a, 
   float *lUn, float *rUn, float *nUn, float *nUn0, float *dx, float *dy, 
   float *dz, float *BCvalsXl, float *BCvalsXr, float *BCvalsYl, float *BCvalsYr, 
   float *BCvalsZl, float *BCvalsZr, int *BCxl, int *BCxr, int *BCyl, int *BCyr, 
   int *BCzl, int *BCzr, int *x0s, int *x0e, int *x1s, int *x1e, int *x2s, 
   int *x2e, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *xlface, int *xrface, 
   int *ylface, int *yrface, int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(mfprob_radresid)(
   float *res1, float *res2, float *res3, float *E1, float *E1old, float *E2, 
   float *E2old, float *E3, float *E3old, float *HI, float *HIold, float *HeI, 
   float *HeIold, float *HeII, float *HeIIold, float *E1src, float *E2src, 
   float *E3src, int *LimType, int *LimImp, float *dt, float *theta, FLOAT *a, 
   FLOAT *a0, FLOAT *adot, FLOAT *adot0, float *aUn, float *lUn, float *lUn0, 
   float *rUn, float *rUn0, float *nUn, float *nUn0, float *dx, float *dy, 
   float *dz, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);
  
extern "C" void FORTRAN_NAME(mfprob_locresid)(
   float *ecres, float *HIres, float *HeIres, float *HeIIres, float *ecsrc, 
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *vx, float *vy, float *vz, 
   float *rho, float *ec, float *eh, float *HI, float *HI0, float *HeI, float *HeI0, 
   float *HeII, float *HeII0, float *piHI, float *piHeI, float *piHeII, float *GHI, 
   float *GHeI, float *GHeII, float *dt, float *theta, int *DualEnergy, FLOAT *a, 
   FLOAT *a0, FLOAT *adot, FLOAT *adot0, float *gamma, float *hfrac, int *model, 
   float *CompA, float *Comp_xray, float *Comp_temp, int *NTempBins, float *TempStart, 
   float *TempEnd, float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, 
   float *k5Tb, float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, 
   float *ciHITb, float *ciHeITb, float *ciHeISTb, float *ciHeIITb, 
   float *reHIITb, float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, 
   float *bremTb, float *aUn, float *dUn, float *dUn0, float *vUn, float *lUn, 
   float *lUn0, float *rUn, float *rUn0, float *eUn, float *nUn, float *nUn0, 
   float *eScale, int *Nchem, float *dx, float *dy, float *dz, int *Nx, int *Ny, 
   int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfprob_radjac)(
   float *E1jac_E1, float *E1jac_HI, float *E2jac_E2, float *E2jac_HI, 
   float *E2jac_HeI, float *E3jac_E3, float *E3jac_HI, float *E3jac_HeI, 
   float *E3jac_HeII, float *E1, float *E2, float *E3, float *HI, float *HeI, 
   float *HeII, float *dt, float *theta, FLOAT *a, FLOAT *adot, float *aUn, 
   float *lUn, float *rUn, float *nUn, int *Nchem, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(blocksolve)(
   float *Amat, float *xvec, float *bvec, int *N, int *M, int *ier);

extern "C" void FORTRAN_NAME(mfprob_sources)(
   float *eta1, float *eta2, float *eta3, float *FS_NGammaDot, float *ecsrc, 
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *time, FLOAT *a, 
   int *Model, int *ProblemType, int *ESpectrum, int *Nchem, int *SMEmiss, 
   float *NGammaDot, float *EtaRadius, float *EtaCenter, float *aUn, 
   float *dUn, float *vUn, float *lUn, float *tUn, float *rUn, float *eUn, 
   float *nUn, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, 
   float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(mfprob_radinit)(
   float *Ef, float *E1, float *E2, float *E3, float *EiScale, 
   int *ESpectrum, int *ier);

extern "C" void FORTRAN_NAME(mfprob_computeradiationintegrals)(
   float *piHI, float *piHeI, float *piHeII, float *GHI, float *GHeI, 
   float *GHeII, float *Ef, float *E1, float *E2, float *E3, float *E1old, 
   float *E2old, float *E3old, int *Nchem, int *ESpectrum, float *fsUn, 
   float *rUn, float *rUn0, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfprob_enforceradiationbounds)(
   float *Ef, float *E1, float *E2, float *E3, float *fsUn, float *rUn, 
   int *ESpectrum, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfprob_analyticresid)(
   float *ecres, float *HIres, float *HeIres, float *HeIIres, float *ec, 
   float *HI, float *HeI, float *HeII, float *HIold, float *HeIold, float *HeIIold, 
   float *dt, float *vx, float *vy, float *vz, float *rho, float *eh, float *ecsrc, 
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *piHI, float *piHeI, 
   float *piHeII, float *GHI, float *GHeI, float *GHeII, float *gamma, float *HFrac, 
   int *model, int *DualEnergy, FLOAT *a, FLOAT *a0, FLOAT *adot, FLOAT *adot0, 
   float *CompA, float *Comp_xray, float *Comp_temp, int *NTempBins, float *TStart, 
   float *TEnd, float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, 
   float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, 
   float *ciHeITb, float *ciHeISTb, float *ciHeIITb, float *reHIITb, 
   float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, float *bremTb, float *aUn, 
   float *dUn, float *dUn0, float *vUn, float *lUn, float *lUn0, float *rUn, 
   float *rUn0, float *eUn, float *nUn, float *nUn0, float *ecScale, int *Nchem, 
   int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(mfprob_analyticinitguess)(
   float *Ef, float *E1, float *E2, float *E3, float *ec, float *HI, float *HeI, 
   float *HeII, float *dt, float *vx, float *vy, float *vz, float *rho, float *eh, 
   float *E1src, float *E2src, float *E3src, float *ecsrc, float *HIsrc, 
   float *HeIsrc, float *HeIIsrc, float *gamma, float *hfrac, int *model, 
   int *dualenergy, int *ESpectrum, FLOAT *a, FLOAT *a0, FLOAT *adot, FLOAT *adot0, 
   float *CompA, float *Comp_xray, float *Comp_temp, float *piHI, float *piHeI, 
   float *piHeII, float *GHI, float *GHeI, float *GHeII, int *NTempBins, float *TStart, 
   float *TEnd, float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, 
   float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, 
   float *ciHeITb, float *ciHeISTb, float *ciHeIITb, float *reHIITb, 
   float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, float *bremTb, float *aUn, 
   float *dUn, float *dUn0, float *vUn, float *lUn, float *lUn0, float *fsUn, 
   float *rUn, float *rUn0, float *eUn, float *nUn, float *nUn0, float *ecScale, 
   int *Nchem, float *dx, float *dy, float *dz, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);




/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the MFProb class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int MFProb::SetupSystem(Eflt64 *mat11, Eflt64 *mat12, Eflt64 *mat13, Eflt64 *mat21, 
			Eflt64 *mat22, Eflt64 *mat23, Eflt64 *mat31, Eflt64 *mat32, 
			Eflt64 *mat33, Eflt64 *rhs1, Eflt64 *rhs2, Eflt64 *rhs3, 
			EnzoVector *u, float *adj11, float *adj12, float *adj13, 
			float *adj21, float *adj22, float *adj23, 
			float *adj31, float *adj32, float *adj33) 
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
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *E1old = U0->GetData(iE1);
  float *E2old = U0->GetData(iE2);
  float *E3old = U0->GetData(iE3);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  FORTRAN_NAME(mfprob_setupsystem)
    (mat11, mat12, mat13, mat21, mat22, mat23, mat31, mat32, mat33, rhs1, rhs2, rhs3, 
     E1, E2, E3, E1old, E2old, E3old, HI, HIold, HeI, HeIold, HeII, HeIIold, adj11, 
     adj12, adj13, adj21, adj22, adj23, adj31, adj32, adj33, &LimType, &LimImp, &dt, 
     &theta, &a, &LenUnits, &rUnits, &nUnits, &nUnits0, &dx[0], &dx[1], &dx[2], 
     EBdryVals[0][0], EBdryVals[0][1], EBdryVals[1][0], EBdryVals[1][1], 
     EBdryVals[2][0], EBdryVals[2][1], &(BdryType[0][0]), &(BdryType[0][1]), 
     &(BdryType[1][0]), &(BdryType[1][1]), &(BdryType[2][0]), &(BdryType[2][1]), &x0s, 
     &x0e, &x1s, &x1e, &x2s, &x2e, &Nchem, &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &xlface, &xrface, &ylface, &yrface, &zlface, &zrface, &ier);
  float ftime = MPI_Wtime();
  timers[0] += ftime - stime;
  return(ier);
}

/********/
int MFProb::RadResid(EnzoVector *radresid, EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *res1 = radresid->GetData(iE1);
  float *res2 = radresid->GetData(iE2);
  float *res3 = radresid->GetData(iE3);
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *E1old = U0->GetData(iE1);
  float *E2old = U0->GetData(iE2);
  float *E3old = U0->GetData(iE3);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  float *E1src = extsrc->GetData(iE1);
  float *E2src = extsrc->GetData(iE2);
  float *E3src = extsrc->GetData(iE3);
  FORTRAN_NAME(mfprob_radresid)
    (res1, res2, res3, E1, E1old, E2, E2old, E3, E3old, HI, HIold, HeI, 
     HeIold, HeII, HeIIold, E1src, E2src, E3src, &LimType, &LimImp, &dt, 
     &theta, &a, &a0, &adot, &adot0, &aUnits, &LenUnits, &LenUnits0, 
     &rUnits, &rUnits0, &nUnits, &nUnits0, &dx[0], &dx[1], &dx[2], &Nchem,
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[1] += ftime - stime;
  return(ier);
}


/********/
int MFProb::LocResid(EnzoVector *locresid, EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *ecres = locresid->GetData(iec);
  float *HIres = locresid->GetData(iHI);
  float *HeIres = locresid->GetData(iHeI);
  float *HeIIres = locresid->GetData(iHeII);
  float *ecsrc = extsrc->GetData(iec);
  float *HIsrc = extsrc->GetData(iHI);
  float *HeIsrc = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  float *ec = u->GetData(iec);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  FORTRAN_NAME(mfprob_locresid)
    (ecres, HIres, HeIres, HeIIres, ecsrc, HIsrc, HeIsrc, HeIIsrc, vx, vy, vz, 
     rho, ec, eh, HI, HIold, HeI, HeIold, HeII, HeIIold, piHI, piHeI, piHeII, 
     GHI, GHeI, GHeII, &dt, &theta, &dualenergy, &a, &a0, &adot, &adot0, &Gamma, 
     &HFrac, &Model, &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &aUnits, &DenUnits, &DenUnits0, 
     &VelUnits, &LenUnits, &LenUnits0, &rUnits, &rUnits0, &eUnits, &nUnits, 
     &nUnits0, &eScale, &Nchem, &dx[0], &dx[1], &dx[2], &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[2] += ftime - stime;
  return(ier);
}


/********/
int MFProb::RadJac(EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  FORTRAN_NAME(mfprob_radjac)
    (L_E1_E1, L_E1_HI, L_E2_E2, L_E2_HI, L_E2_HeI, L_E3_E3, L_E3_HI, 
     L_E3_HeI, L_E3_HeII, E1, E2, E3, HI, HeI, HeII, &dt, &theta, &a, &adot, 
     &aUnits, &LenUnits, &rUnits, &nUnits, &Nchem, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[3] += ftime - stime;
  return(ier);
}


/********/
int MFProb::BlockSolve(float *Amat, float *xvec, float *bvec, int *N, int *M)
{
  float stime = MPI_Wtime();
  int ier;
  FORTRAN_NAME(blocksolve)(Amat, xvec, bvec, N, M, &ier);
  float ftime = MPI_Wtime();
  timers[4] += ftime - stime;
  return(ier);
}


/********/
int MFProb::Sources(EnzoVector *extsrc, float *time, float *FS_NGammaDot)
{
  float stime = MPI_Wtime();
  int ier;
  float *eta1 = extsrc->GetData(iE1);
  float *eta2 = extsrc->GetData(iE2);
  float *eta3 = extsrc->GetData(iE3);
  float *ecsrc = extsrc->GetData(iec);
  float *HIsrc = extsrc->GetData(iHI);
  float *HeIsrc = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  FORTRAN_NAME(mfprob_sources)
    (eta1, eta2, eta3, FS_NGammaDot, ecsrc, HIsrc, HeIsrc, HeIIsrc, time, 
     &a, &Model, &ProblemType, &ESpectrum, &Nchem, &StarMakerEmissivityField, 
     &IonizationParms[0], &IonizationParms[1], &IonizationParms[2], 
     &aUnits, &DenUnits, &VelUnits, &LenUnits, &TimeUnits, &rUnits, 
     &eUnits, &nUnits, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[5] += ftime - stime;
  return(ier);
}


/********/
int MFProb::RadInit(float *Ef, float *E1, float *E2, float *E3, float *EiScale)
{
  int ier;
  FORTRAN_NAME(mfprob_radinit)(Ef, E1, E2, E3, EiScale, &ESpectrum, &ier);
  return(ier);
}


/********/
int MFProb::ComputeRadiationIntegrals(EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *E1 = u->GetData(iE1);
  float *E2 = u->GetData(iE2);
  float *E3 = u->GetData(iE3);
  float *E1old = U0->GetData(iE1);
  float *E2old = U0->GetData(iE2);
  float *E3old = U0->GetData(iE3);
  FORTRAN_NAME(mfprob_computeradiationintegrals)
    (piHI, piHeI, piHeII, GHI, GHeI, GHeII, Efree, E1, E2, E3, E1old, E2old, 
     E3old, &Nchem, &ESpectrum, &fsUnits, &rUnits, &rUnits0, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[6] += ftime - stime;
  return(ier);
}


/********/
int MFProb::EnforceRadiationBounds(float *Ef, float *E1, float *E2, float *E3)
{
  float stime = MPI_Wtime();
  int ier;
  FORTRAN_NAME(mfprob_enforceradiationbounds)
    (Ef, E1, E2, E3, &fsUnits, &rUnits, &ESpectrum, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[7] += ftime - stime;
  return(ier);
}


/********/
int MFProb::AnalyticResid(EnzoVector *fu, EnzoVector *u)
{
  float stime = MPI_Wtime();
  int ier;
  float *ecres = fu->GetData(iec);
  float *HIres = fu->GetData(iHI);
  float *HeIres = fu->GetData(iHeI);
  float *HeIIres = fu->GetData(iHeII);
  float *ec = u->GetData(iec);
  float *HI = u->GetData(iHI);
  float *HeI = u->GetData(iHeI);
  float *HeII = u->GetData(iHeII);
  float *HIold = U0->GetData(iHI);
  float *HeIold = U0->GetData(iHeI);
  float *HeIIold = U0->GetData(iHeII);
  float *ecsrc = extsrc->GetData(iec);
  float *HIsrc = extsrc->GetData(iHI);
  float *HeIsrc = extsrc->GetData(iHeI);
  float *HeIIsrc = extsrc->GetData(iHeII);
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  FORTRAN_NAME(mfprob_analyticresid)
    (ecres, HIres, HeIres, HeIIres, ec, HI, HeI, HeII, HIold, HeIold, HeIIold, 
     &dt, vx, vy, vz, rho, eh, ecsrc, HIsrc, HeIsrc, HeIIsrc, piHI, piHeI, piHeII, 
     GHI, GHeI, GHeII, &Gamma, &HFrac, &Model, &dualenergy, &a, &a0, &adot, &adot0, 
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &aUnits, &DenUnits, &DenUnits0, 
     &VelUnits, &LenUnits, &LenUnits0, &rUnits, &rUnits0, &eUnits, &nUnits, 
     &nUnits0, &eScale, &Nchem, &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[8] += ftime - stime;
  return(ier);
}


/********/
int MFProb::AnalyticInitGuess(EnzoVector *u)
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
  FORTRAN_NAME(mfprob_analyticinitguess)
    (Efree, E1, E2, E3, ec, HI, HeI, HeII, &dt, vx, vy, vz, rho, eh, 
     E1src, E2src, E3src, ecsrc, HIsrc, HeIsrc, HeIIsrc, &Gamma, 
     &HFrac, &Model, &dualenergy, &ESpectrum, &a, &a0, &adot, &adot0, 
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, piHI, piHeI, 
     piHeII, GHI, GHeI, GHeII, &CoolData.NumberOfTemperatureBins, 
     &CoolData.TemperatureStart, &CoolData.TemperatureEnd, RateData.k1, 
     RateData.k2, RateData.k3, RateData.k4, RateData.k5, RateData.k6, 
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI, 
     CoolData.ciHeI, CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
     CoolData.reHeII1, CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &aUnits, 
     &DenUnits, &DenUnits0, &VelUnits, &LenUnits, &LenUnits0, &fsUnits, &rUnits, 
     &rUnits0, &eUnits, &nUnits, &nUnits0, &eScale, &Nchem, &dx[0], &dx[1], 
     &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  float ftime = MPI_Wtime();
  timers[9] += ftime - stime;
  return(ier);
}


#endif
