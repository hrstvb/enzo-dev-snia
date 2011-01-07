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
/  Split Implicit Problem Class, Destructor routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



AMRFLDSplit::~AMRFLDSplit()
{

//   if (debug)  printf("Entering AMRFLDSplit::destructor routine\n");

  // delete boundary condition arrays
  int i, j;
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) 
      if (BdryVals[i][j] != NULL)  
	delete [] BdryVals[i][j];


}
#endif
