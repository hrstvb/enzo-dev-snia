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
/  Split Implicit Problem Class, FillRates routine.
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: Fills the photo-ionization and photo-heating arrays used by
/           chemistry and cooling routines using time-averaged internal
/           values for the radiation.
/ 
/  NOTE: In order to save on memory, Enzo's chemistry routines assume 
/        that the photo-heating rates are combined into a single rate, 
/        and scaled by the current number density of HI, to be later 
/        unpacked by rescaling back with HI.  This potentially loses 
/        accuracy in the case that during chemistry subcycling the 
/        HI, HeI and HeII values change significantly, since we retain 
/        the initial rate scaling but use updated HI values in 
/        chemistry subcycles.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"


int AMRFLDSplit::FillRates(LevelHierarchyEntry *LevelArray[], int level)
{

//   if (debug)
//     printf("Entering AMRFLDSplit::FillRates routine\n");

  // set some physical constants
  float c = 2.99792458e10;        // speed of light [cm/s]
  float hp = 6.6260693e-27;       // Planck's constant [ergs*s]
  float mp = 1.67262171e-24;      // mass of a proton [g]
  float ev2erg = 1.60217653e-12;  // conversion constant from eV to ergs
  float dom = DenUnits*POW(a,3.0)/mp;
  float tbase1 = TimeUnits;
  float xbase1 = LenUnits/a/aUnits;
  float dbase1 = DenUnits*POW(a*aUnits,3.0);
  float coolunit = POW(aUnits,5.0) * POW(xbase1,2.0) * POW(mp,2.0) 
    / POW(tbase1,3.0) / dbase1;
  float rtunits = ev2erg/TimeUnits/coolunit/dom;
  
  // iterate over grids owned by this processor (this level down)
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	// set grid dimension information
	int dim, i;
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	          - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;
	int size = x0len*x1len*x2len;
	

	// access Enzo fields 
	float *Enew       = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	float *phHI       = Temp->GridHierarchyEntry->GridData->AccessKPhHI();
	float *phHeI      = Temp->GridHierarchyEntry->GridData->AccessKPhHeI();
	float *phHeII     = Temp->GridHierarchyEntry->GridData->AccessKPhHeII();
	float *photogamma = Temp->GridHierarchyEntry->GridData->AccessPhotoGamma();
	float *dissH2I    = Temp->GridHierarchyEntry->GridData->AccessKDissH2I();
	float *HI         = Temp->GridHierarchyEntry->GridData->AccessHIDensity();
	float *HeI=NULL, *HeII=NULL;
	if (Nchem == 3) {
	  HeI  = Temp->GridHierarchyEntry->GridData->AccessHeIDensity();
	  HeII = Temp->GridHierarchyEntry->GridData->AccessHeIIDensity();
	}

	// check that field data exists
	if (Enew == NULL)
	  ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");
	if (phHI == NULL)
	  ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");
	if (photogamma == NULL)
	  ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");
	if (RadiativeTransferHydrogenOnly == FALSE) {
	  if (phHeI == NULL)
	    ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");
	  if (phHeII == NULL)
	    ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");
	}
	if (MultiSpecies > 1) 
	  if (dissH2I == NULL)
	    ENZO_FAIL("AMRFLDSplit_FillRates ERROR: no radiation array!");

	// fill HI photo-ionization rate
	float pHIconst = c*TimeUnits*intSigESigHInu/hp/intSigE;
	for (i=0; i<size; i++)  phHI[i] = Enew[i]*ErUnits*pHIconst;

	// fill HeI and HeII photo-ionization rates
	float pHeIconst  = c*TimeUnits*intSigESigHeInu/hp/intSigE;
	float pHeIIconst = c*TimeUnits*intSigESigHeIInu/hp/intSigE;
	if (RadiativeTransferHydrogenOnly == FALSE) {
	  for (i=0; i<size; i++)  phHeI[i]  = Enew[i]*ErUnits*pHeIconst;
	  for (i=0; i<size; i++)  phHeII[i] = Enew[i]*ErUnits*pHeIIconst;
	}
   
	// fill photo-heating rate
	float phScale    = c*TimeUnits/intSigE/VelUnits/VelUnits/mp/rtunits;
	float GHIconst   = phScale*(intSigESigHI   - 13.6*ev2erg/hp*intSigESigHInu);
	float GHeIconst  = phScale*(intSigESigHeI  - 24.6*ev2erg/hp*intSigESigHeInu);
	float GHeIIconst = phScale*(intSigESigHeII - 54.4*ev2erg/hp*intSigESigHeIInu);
	if (Nchem == 1)
	  for (i=0; i<size; i++)  photogamma[i] = Enew[i]*ErUnits*GHIconst;
	if (Nchem == 3)
	  for (i=0; i<size; i++)  
	    photogamma[i] = Enew[i]*ErUnits *
	      (GHIconst*HI[i] + GHeIconst*HeI[i] + GHeIIconst*HeII[i])/HI[i];

	// fill H2 dissociation rate (none for grey FLD problems)
	if (MultiSpecies > 1) 
	  for (i=0; i<size; i++)  dissH2I[i] = 0.0;

      }  // end iteration over grids on this processor

  // return success
  return SUCCESS;

}
#endif
