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
/  Destructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit MF problem.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"



MFProb::~MFProb()
{

//   if (debug)  printf("Entering MFProb::destructor routine\n");

  // delete free-streaming solver objects
  delete FSSolve;

  // delete HYPRE objects
#ifdef USE_HYPRE
  HYPRE_SStructStencilDestroy(stencil1);
  HYPRE_SStructStencilDestroy(stencil2);
  HYPRE_SStructStencilDestroy(stencil3);
  HYPRE_SStructGridDestroy(grid);
  HYPRE_SStructGraphDestroy(graph);
  HYPRE_SStructVectorDestroy(rhsvec);
  HYPRE_SStructVectorDestroy(solvec);
  HYPRE_SStructMatrixDestroy(P);
#endif
  if (P11tmpvec != NULL)  delete[] P11tmpvec;
  if (P12tmpvec != NULL)  delete[] P12tmpvec;
  if (P13tmpvec != NULL)  delete[] P13tmpvec;
  if (P21tmpvec != NULL)  delete[] P21tmpvec;
  if (P22tmpvec != NULL)  delete[] P22tmpvec;
  if (P23tmpvec != NULL)  delete[] P23tmpvec;
  if (P31tmpvec != NULL)  delete[] P31tmpvec;
  if (P32tmpvec != NULL)  delete[] P32tmpvec;
  if (P33tmpvec != NULL)  delete[] P33tmpvec;
  if (r1tmpvec  != NULL)  delete[] r1tmpvec;
  if (r2tmpvec  != NULL)  delete[] r2tmpvec;
  if (r3tmpvec  != NULL)  delete[] r3tmpvec;
  if (HYPREbuff != NULL)  delete[] HYPREbuff;

  // delete InexactNewton solver object
  delete INSolve;

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
  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete L_e;
  delete L_HI;
  delete L_HeI;
  delete L_HeII;

  //   arrays require deleting the array
  if (L_E1_E1   != NULL)  delete[] L_E1_E1;
  if (L_E2_E2   != NULL)  delete[] L_E2_E2;
  if (L_E3_E3   != NULL)  delete[] L_E3_E3;
  if (L_E1_HI   != NULL)  delete[] L_E1_HI;
  if (L_E2_HI   != NULL)  delete[] L_E2_HI;
  if (L_E2_HeI  != NULL)  delete[] L_E2_HeI;
  if (L_E3_HI   != NULL)  delete[] L_E3_HI;
  if (L_E3_HeI  != NULL)  delete[] L_E3_HeI;
  if (L_E3_HeII != NULL)  delete[] L_E3_HeII;
  if (eCorr     != NULL)  delete[] eCorr;

}
#endif
