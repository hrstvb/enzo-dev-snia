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
/  Split Implicit Problem Class, Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Computes the radiation time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
 
 
float AMRFLDSplit::ComputeTimeStep(LevelHierarchyEntry *LevelArray[], int level)
{

  // Set time step depending on how it has been set up by the user:
  //    If dtfac is set, compute maximum time step as estimate 
  //    allowing dtfac relative change.  This relative change is
  //    estimated as follows:
  //       dt_new = dt_old / relerr_fac
  //    where relerr_fac gives the ratio between an estimated 
  //    local truncation error and the desired relative change:
  //       relerr_fac = || (unew - uold) / w ||_p
  //    with the scaling vector w given by
  //       w = dtfac*[sqrt(|unew*uold|) + atol]
  //    and where we have the following parameters:
  //       p - norm choice (input), 0->max norm, otherwise the p-norm
  //           **all p-norms here divide by the total number of cells**
  //       dtfac - desired relative change per step (input)
  //       atol - 0.1 (assumes units are all normalized)
  float dt_est = huge_number;    // max time step (normalized units)
  if (dtfac != huge_number) {

    // initialize variables
    float diff, w, tmp, atol;
    int i, j, k, l;
    float loc_est;

    // iterate over grids owned by this processor (this level down)
    LevelHierarchyEntry* Temp;
    for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
      for (Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
	if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
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

	  // set a cell "volume" assuming the global domain has volume 1
	  float dV = 1;
	  for (int dim=0; dim<rank; dim++)
	    dV *= (Temp->GridData->GetGridRightEdge(dim) 
		   - Temp->GridData->GetGridLeftEdge(dim)) 
	        / (Temp->GridData->GetGridDimension(dim))
	        / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
      
      
	  // access old/new radiation fields (old stored in KPhHI)
	  float *Eold = Temp->GridData->AccessKPhHI();
	  float *Enew = Temp->GridData->AccessRadiationFrequency0();
	  
	  // perform estimates for the radiation energy relative change on this grid
	  atol = 0.1;   // assumes values are normalized
	  if (dtnorm > 0.0) {
	    for (k=ghZl; k<n3[2]+ghZl; k++) 
	      for (j=ghYl; j<n3[1]+ghYl; j++)
		for (i=ghXl; i<n3[0]+ghXl; i++) {
		  w = dtfac*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				       *Eold[(k*x1len + j)*x0len + i])) 
			     + atol);
		  diff = Enew[(k*x1len + j)*x0len + i] 
		       - Eold[(k*x1len + j)*x0len + i];
		  tmp = fabs(diff/w);
		  loc_est += POW(tmp,dtnorm)/dV;
		}
	  }
	  else {
	    for (k=ghZl; k<n3[2]+ghZl; k++) 
	      for (j=ghYl; j<n3[1]+ghYl; j++)
		for (i=ghXl; i<n3[0]+ghXl; i++) {
		  w = dtfac*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				       *Eold[(k*x1len + j)*x0len + i])) 
			     + atol);
		  diff = Enew[(k*x1len + j)*x0len + i] 
		       - Eold[(k*x1len + j)*x0len + i];
		  tmp = fabs(diff/w);
		  loc_est = (loc_est > tmp) ? loc_est : tmp;
		}
	  }
	  
	}  // end iteration over grids on this processor

    // communicate to obtain overall sum/max
    float glob_est;
#ifdef USE_MPI
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg ONE = 1;
    if (dtnorm > 0.0) 
      MPI_Allreduce(&loc_est,&glob_est,ONE,DataType,MPI_SUM,MPI_COMM_WORLD);
    else
      MPI_Allreduce(&loc_est,&glob_est,ONE,DataType,MPI_MAX,MPI_COMM_WORLD);
#else
    glob_est = loc_est;
#endif

    // compute overall norms
    if (dtnorm > 0.0)  glob_est = POW(glob_est,1.0/dtnorm);

    // compute time step estimate (physical units)
    dt_est = (glob_est == 0.0) ? huge_number : dt/glob_est;
    dt_est = min(dt_est, huge_number);

  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;
}

#endif
