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
/  Problem parameter output routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"

int MFProb::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering MFProb::WriteParameters routine\n");
  
  fprintf(fptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
  fprintf(fptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "RadHydroHFraction = %"FSYM"\n", HFrac);
  fprintf(fptr, "RadHydroModel = %"ISYM"\n", Model);

  // set restart initial time step to current time step
  fprintf(fptr, "RadHydroMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "RadHydroMinDt = %22.16e\n", mindt);
  if (dt == 0.0) 
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", initdt);
  else
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", dt);
  fprintf(fptr, "RadHydroDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "RadHydroDtRadFac = %22.16e\n", dtfac[0]);
  fprintf(fptr, "RadHydroDtGasFac = %22.16e\n", dtfac[1]);
  fprintf(fptr, "RadHydroDtChemFac = %22.16e\n", dtfac[2]);
  fprintf(fptr, "RadiationScaling = %22.16e\n", rScale);
  fprintf(fptr, "EnergyCorrectionScaling = %22.16e\n", eScale);
  fprintf(fptr, "ChemistryScaling = %22.16e\n", nScale);
  fprintf(fptr, "RadHydroTheta = %22.16e\n", theta);
  fprintf(fptr, "RadHydroLimiterType = %"ISYM"\n", LimType);
  fprintf(fptr, "RadHydroImplicitLimiter = %"ISYM"\n", LimImp);
  fprintf(fptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) 
    fprintf(fptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[1][0], BdryType[1][1]);
  if (rank > 2) 
    fprintf(fptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[2][0], BdryType[2][1]);
  fprintf(fptr, "RadHydroAprxJacobian = %"ISYM"\n", approx_jac);    
  fprintf(fptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);    
  fprintf(fptr, "RadHydroAnalyticChem = %"ISYM"\n", AnalyticChem);
  fprintf(fptr, "RadHydroSemiImplicit = %"ISYM"\n", semi_implicit);
  fprintf(fptr, "RadHydroNewtLinesearch = %"ISYM"\n", newt_linesearch);    
  fprintf(fptr, "RadHydroNewtIters = %"ISYM"\n", newt_maxit);    
  fprintf(fptr, "RadHydroNewtNorm = %"ISYM"\n", newt_norm);    
  fprintf(fptr, "RadHydroINConst = %22.16e\n", newt_INconst);    
  fprintf(fptr, "RadHydroNewtTolerance = %22.16e\n", newt_tol);    
  fprintf(fptr, "RadHydroMinLinesearch = %22.16e\n", newt_MinLinesearch);    
  fprintf(fptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "RadHydroMGPostRelax = %i\n", sol_npost);    

  // if doing an ionization problem (ProblemTypes 410-415),  
  // output additional parameters 
  if ((ProblemType >= 410) && (ProblemType <= 415)) {
    fprintf(fptr, "NGammaDot = %22.16e\n", IonizationParms[0]);
    fprintf(fptr, "EtaRadius = %22.16e\n", IonizationParms[1]);
    fprintf(fptr, "EtaCenter = %22.16e %22.16e %22.16e\n",
	    IonizationParms[2],IonizationParms[3],IonizationParms[4]);
  }
  
  return SUCCESS;
}
#endif
