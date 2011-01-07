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
/  Split Implicit Problem Class, Dump routine
/  
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Outputs entire problem state to stderr and files.  Useful 
/           upon solve failure.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"


int AMRFLDSplit::Dump(EnzoVector *ucur)
{

  if (debug) {
    fprintf(stderr,"Dumping AMRFLDSplit module parameters to file RTdump.params\n");
    FILE *fptr = fopen("RTdump.params", "w");
    this->WriteParameters(fptr);
    fclose(fptr);
  }

  return SUCCESS;
}
#endif
