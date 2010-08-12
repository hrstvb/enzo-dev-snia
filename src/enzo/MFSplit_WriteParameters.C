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
#include "MFSplit.h"


int MFSplit::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering MFSplit::WriteParameters routine\n");
  
  fprintf(fptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
  fprintf(fptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "RadHydroHFraction = %22.16e\n", HFrac);
  fprintf(fptr, "RadHydroModel = %"ISYM"\n", Model);

  // set restart initial time step to current time step
  if (dt == 0.0) 
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", initdt);
  else
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", dt);
  fprintf(fptr, "RadHydroMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "RadHydroMinDt = %22.16e\n", mindt);
  fprintf(fptr, "RadHydroDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "RadHydroDtRadFac = %22.16e\n", dtfac[0]);
  fprintf(fptr, "RadHydroDtGasFac = %22.16e\n", dtfac[1]);
  fprintf(fptr, "RadHydroDtChemFac = %22.16e\n", dtfac[2]);

  fprintf(fptr, "RadiationScaling = %22.16e\n", rScale);
  fprintf(fptr, "EnergyCorrectionScaling = %22.16e\n", eScale);
  fprintf(fptr, "ChemistryScaling = %22.16e\n", nScale);

  fprintf(fptr, "E1Units = %22.16e\n", E1Units*rUnits);
  fprintf(fptr, "E2Units = %22.16e\n", E2Units*rUnits);
  fprintf(fptr, "E3Units = %22.16e\n", E3Units*rUnits);

  fprintf(fptr, "RadHydroTheta = %22.16e\n", theta);
  fprintf(fptr, "RadHydroLimiterType = %"ISYM"\n", LimType);

  fprintf(fptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) 
    fprintf(fptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[1][0], BdryType[1][1]);
  if (rank > 2) 
    fprintf(fptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[2][0], BdryType[2][1]);

  fprintf(fptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);    
  fprintf(fptr, "RadHydroNewtTolerance = %22.16e\n", sol_tolerance);    
  fprintf(fptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "RadHydroMGPostRelax = %i\n", sol_npost);    

  fprintf(fptr, "NGammaDot = %22.16e\n", NGammaDot);
  fprintf(fptr, "EtaRadius = %22.16e\n", EtaRadius);
  fprintf(fptr, "EtaCenter = %22.16e %22.16e %22.16e\n",
	  EtaCenter[0], EtaCenter[1], EtaCenter[2]);

  return SUCCESS;
}
#endif
