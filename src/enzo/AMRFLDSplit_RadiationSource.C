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
/  Split Implicit Problem Class, Emissivity Field Computation Routine
/
/  written by: Daniel Reynolds
/  date:       January 2011
/  modified1:  
/
/  PURPOSE: Computes the emissivity field to enforce on the radiation 
/           energy equation.  This is only called if StarMaker sources 
/           are not handled elsewhere (i.e. for test problems).
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
 
 
float AMRFLDSplit::RadiationSource(LevelHierarchyEntry *LevelArray[], 
				   int level, float time)
{

//   if (debug)  printf("    AMRFLDSplit RadiationSource, initializing local variables\n");

  // initialize local variables to be reused
  int i, j, k;
  float SpecConst = (ESpectrum == 1) ? 1.52877652583602 : 1.0;
  float h_nu0 = 13.6 * 1.60217653e-12;
  float total_eta = 0.0;
  float etaconst;
  float cellZl, cellZr, cellYl, cellYr, cellXl, cellXr, cellXc, cellYc, cellZc;
  
//   if (debug)  printf("    AMRFLDSplit RadiationSource, iterating over processor grids\n");

  // iterate over grids owned by this processor (this level down)
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber()) {
	
	// set dimension information
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridData->GetGridEndIndex(dim)
  	          - Temp->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;
	float x0L = Temp->GridData->GetGridLeftEdge(0);
	float x1L = Temp->GridData->GetGridLeftEdge(1);
	float x2L = Temp->GridData->GetGridLeftEdge(2);
	float dx[3] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  dx[dim] = (Temp->GridData->GetGridRightEdge(dim) 
		     - Temp->GridData->GetGridLeftEdge(dim)) 
	          / n3[dim];
	float dV = dx[0]*dx[1]*dx[2]*LenUnits*LenUnits*LenUnits;

	// set a cell "normalized volume" assuming the global domain has volume 1
	float dVscale = 1;
	for (int dim=0; dim<rank; dim++)
	  dVscale *= (Temp->GridData->GetGridRightEdge(dim) 
 		    - Temp->GridData->GetGridLeftEdge(dim)) 
	            / n3[dim] / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	
	// access emissivity field, initialize to zero
	float *eta = Temp->GridData->AccessEmissivity0();
	for (i=0; i<x0len*x1len*x2len; i++)  eta[i] = 0.0;
	
	
	/////////////
	// compute emissivity field based on ProblemType, input parameters
	
	// source with a fixed location/radius
	if ( (ProblemType == 410) || (ProblemType == 411) || 
	     (ProblemType == 413) || (ProblemType == 415) ) {
	  
	  // one-cell source
	  if (EtaRadius == 0.0) {
	    
	    // compute eta factor for given ionization source
	    etaconst = h_nu0 * NGammaDot * SpecConst / dV;

	    for (k=ghZl; k<n3[2]+ghZl; k++) {
	      
	      // z-boundaries (comoving) for this cell
	      cellZl = x2L + (k-ghZl)*dx[2];
	      cellZr = cellZl + dx[2];
	      
	      // check if we're in the right z-cell
	      if ((cellZl > EtaCenter[2]) || (cellZr <= EtaCenter[2]))
		continue; 
	      
	      for (j=ghYl; j<n3[1]+ghYl; j++) {
		
		// y-boundaries (comoving) for this cell
		cellYl = x1L + (j-ghYl)*dx[1];
		cellYr = cellYl + dx[1];
		
		// check if we're in the right y-cell
		if ((cellYl > EtaCenter[1]) || (cellYr <= EtaCenter[1]))
		  continue; 
		
		for (i=ghXl; i<n3[0]+ghXl; i++) {
		  
		  // x-boundaries (comoving) for this cell
		  cellXl = x0L + (i-ghXl)*dx[0];
		  cellXr = cellXl + dx[0];
		  
		  // check if we're in the right x-cell
		  if ( (cellXl <= EtaCenter[0]) && (cellXr > EtaCenter[0])) { 
		    eta[(k*x1len + j)*x0len + i] = etaconst;
		    break;
		  }
		} // x-loop
	      } // y-loop
	    } // z-loop
	    
	  } else { // multi-celled source
	    
	    // compute eta factor for given ionization source
	    etaconst = h_nu0 * NGammaDot * SpecConst / dV / POW(2.0*EtaRadius, rank);
	    
	    for (k=ghZl; k<n3[2]+ghZl; k++) {
	      
	      // z-center (comoving) for this cell
	      cellZc = x2L + (k-ghZl+0.5)*dx[2];
	      
	      for (j=ghYl; j<n3[1]+ghYl; j++) {
		
		// y-center (comoving) for this cell
		cellYc = x1L + (j-ghYl+0.5)*dx[1];
		
		for (i=ghXl; i<n3[0]+ghXl; i++) {
		  
		  // x-center (comoving) for this cell
		  cellXc = x0L + (i-ghXl+0.5)*dx[0];
		  
		  // see if cell is within source region
		  if ( (fabs(cellXc-EtaCenter[0]) < EtaRadius*dx[0]) &&
		       (fabs(cellYc-EtaCenter[1]) < EtaRadius*dx[1]) &&
		       (fabs(cellZc-EtaCenter[2]) < EtaRadius*dx[2]) )
		    eta[(k*x1len + j)*x0len + i] = etaconst;
		  
		} // x-loop
	      } // y-loop
	    } // z-loop
	  } // EtaRadius == 0
	} // ProblemType in {410, 411, 413, 415}
	
	
	// emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
	if (ProblemType == 412) {
	  
	  // place ionization source along left wall (if on this subdomain)
	  if (x0L == 0.0) {
	    
	    // compute eta factor for given ionization source, and put on wall
	    etaconst = h_nu0 * NGammaDot * SpecConst / dx[0];
	    
	    // iterate over wall
	    for (k=ghZl; k<n3[2]+ghZl; k++)
	      for (j=ghYl; j<n3[1]+ghYl; j++)
		eta[(k*x1len + j)*x0len + ghXl] = etaconst;
	  } // if x0L
	} // ProblemType == 412
	
	
	  //  point-source emissivity at center of every processor
	if (ProblemType == 414) {
	  
	  // compute eta factor for given ionization source
	  etaconst = h_nu0 * NGammaDot * SpecConst / dV;
	  
	  // place ionization source in center of subdomain
	  i = x0len / 2;
	  j = x1len / 2;
	  k = x2len / 2;
	  eta[(k*x1len + j)*x0len + ghXl] = etaconst;
	  
	} // ProblemType == 414
	
	
	// homogeneous emissivity field w/ strength hnu0*NGammaDot/dV
	if (ProblemType == 416) {
	  
	  // compute eta factor for given ionization source
	  etaconst = h_nu0 * NGammaDot * SpecConst / dV;
	  for (i=0; i<x0len*x1len*x2len; i++)  eta[i] = etaconst;
	  
	} // ProblemType == 416
	
	
	// L2 norm of total active emissivity on this grid
	for (k=ghZl; k<n3[2]+ghZl; k++) 
	  for (j=ghYl; j<n3[1]+ghYl; j++) 
	    for (i=ghXl; i<n3[0]+ghXl; i++) 
 	      total_eta += eta[(k*x1len+j)*x0len+i]*eta[(k*x1len+j)*x0len+i]*dVscale;

      }  // end iteration over grids on this processor
  
  // communicate to obtain overall emissivity
  float glob_eta;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg ONE = 1;
  MPI_Allreduce(&total_eta,&glob_eta,ONE,DataType,MPI_SUM,MPI_COMM_WORLD);
#else
  glob_eta = total_eta;
#endif
  
  return sqrt(glob_eta);

}

#endif
