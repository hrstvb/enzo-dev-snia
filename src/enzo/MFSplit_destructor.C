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
/  Destructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the split MF problem.
/
************************************************************************/
#ifdef TRANSFER
#include "MFSplit.h"



MFSplit::~MFSplit()
{

//   if (debug)  printf("Entering MFSplit::destructor routine\n");

  // delete free-streaming solver objects
  delete FSSolve;

  // delete HYPRE objects
#ifdef USE_HYPRE
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
#endif
  if (matentries != NULL)  delete[] matentries;
  if (rhsentries != NULL)  delete[] rhsentries;
  if (HYPREbuff  != NULL)  delete[] HYPREbuff;

  // delete boundary condition arrays
  int i, j;
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) 
      if (EBdryVals[i][j] != NULL)  delete[] EBdryVals[i][j];


  // delete EnzoVectors and other internal arrays
  //   EnzoVectors require deleting the structure
  delete sol;
  delete U0;
  delete extsrc;

  //   arrays require deleting the array
  if (eCorr != NULL)  delete[] eCorr;

}
#endif
