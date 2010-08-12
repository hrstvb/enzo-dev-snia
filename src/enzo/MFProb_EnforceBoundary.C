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
/  EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Enforces boundary conditions on a MFProb vector.  
/           Depending on 'flag' this routine will perform
/           one of two BC-related actions:
/
/              flag=0: sets Dirichlet, Neumann, or mixed (Robin) values 
/              into a given vector.  Useful for enforcing conditions on 
/              a solution.
/
/              flag!=0: set zero-valued homogeneous Dirichlet, 
/              Neumann or mixed (Robin) conditions on any external face.  
/              Useful for enforcing conditions into a Newton update, 
/              which should not interfere with BC values.
/
/           Note: Neumann values are enforced on the first layer of 
/                 ghost zones using a first-order central difference 
/                 approximation to the outward-normal derivative.
/           Note: Since the internal radiation variables are comoving 
/                 and normalized, we renormalize the boundary conditions 
/                 as they are enforced to match the internal units.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"



int MFProb::EnforceBoundary(EnzoVector *u, int flag)
{

  float stime = MPI_Wtime();
  
//   if (debug)
//     printf("Entering MFProb::EnforceBoundary routine\n");

  // check that the MFProb has been prepared
  if (!prepared) 
    ENZO_FAIL("MFProb EnforceBoundary: MFProb unprepared\n");
  
  // get information about the vector u, and check against BC dims
  int i, i2, j, j2, k, k2, idx, idx2, idxbc;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  if (udims[0] != LocDims[0]) 
    ENZO_VFAIL("p%"ISYM" EnforceBC: mismatched x0 dims %"ISYM"!=%"ISYM"\n",
	       MyProcessorNumber,udims[0],LocDims[0])
  if (udims[1] != LocDims[1]) 
    ENZO_VFAIL("p%"ISYM" EnforceBC: mismatched x1 dims %"ISYM"!=%"ISYM"\n",
	       MyProcessorNumber,udims[1],LocDims[1])
  if (udims[2] != LocDims[2]) 
    ENZO_VFAIL("p%"ISYM" EnforceBC: mismatched x2 dims %"ISYM"!=%"ISYM"\n",
	       MyProcessorNumber,udims[2],LocDims[2])
  if (udims[3] != (4+Nchem)) 
    ENZO_VFAIL("p%"ISYM" EnforceBC: mismatched nspecies %"ISYM"!=%"ISYM"\n",
	       MyProcessorNumber,udims[3],4+Nchem)

  // set some shortcuts for the EnzoVector dimensions, scaling
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];
  float dxa = dx[0]*LenUnits/a;
  float dya = dx[1]*LenUnits/a;
  float dza = dx[2]*LenUnits/a;
  float *udata;

  // flag == 0: set BCs into udata
  if (flag == 0) {

    //////////////////////////////////////////////
    // first handle radiation energy frequency 1
    udata = u->GetData(iE1);

    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<ugh[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][0][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][0] && (BdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][0][idxbc]/rUnits;
	}
    }

    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (BdryType[0][1]==1))
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][1][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][1] && (BdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][1][idxbc]/rUnits;
	}
    }

    if (rank > 1) {
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][0] && (BdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][0][idxbc]/rUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][1] && (BdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][0] && (BdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][0][idxbc]/rUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][1] && (BdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 2
      


    //////////////////////////////////////////////
    // handle radiation energy frequency 2
    udata = u->GetData(iE2);

    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<ugh[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][0][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][0] && (BdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][0][idxbc]/rUnits;
	}
    }

    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (BdryType[0][1]==1))
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][1][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][1] && (BdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][1][idxbc]/rUnits;
	}
    }

    if (rank > 1) {
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][0] && (BdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][0][idxbc]/rUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][1] && (BdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][0] && (BdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][0][idxbc]/rUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][1] && (BdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 2
      


    //////////////////////////////////////////////
    // handle radiation energy frequency 3
    udata = u->GetData(iE3);

    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<ugh[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][0][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][0] && (BdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][0][idxbc]/rUnits;
	}
    }

    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (BdryType[0][1]==1)) 
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = EBdryVals[0][1][idxbc]/rUnits;
	  }
    //   Neumann
    if (OnBdry[0][1] && (BdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*EBdryVals[0][1][idxbc]/rUnits;
	}
    }

    if (rank > 1) {
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][0] && (BdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][0][idxbc]/rUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = EBdryVals[1][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[1][1] && (BdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*EBdryVals[1][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][0][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][0] && (BdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][0][idxbc]/rUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = EBdryVals[2][1][idxbc]/rUnits;
	    }
      //   Neumann
      if (OnBdry[2][1] && (BdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*EBdryVals[2][1][idxbc]/rUnits;
	  }
      }
    }  // end if rank > 2
      

  }
  // flag != 0: enforce zero values at Dirichlet faces
  else {

    //////////////////////////////////////////////
    // first handle radiation energy frequency 1
    udata = u->GetData(iE1);

    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=0; i<ugh[0][0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0][1]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    if (rank > 1) {
      // x1 left boundary
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x1 right boundary
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    if (rank > 2) {
      // x2 left boundary
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x2 right boundary
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    //////////////////////////////////////////////
    // handle radiation energy frequency 2
    udata = u->GetData(iE2);

    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=0; i<ugh[0][0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0][1]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    if (rank > 1) {
      // x1 left boundary
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x1 right boundary
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    if (rank > 2) {
      // x2 left boundary
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x2 right boundary
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    //////////////////////////////////////////////
    // first handle radiation energy frequency 3
    udata = u->GetData(iE3);

    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0][0]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=0; i<ugh[0][0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0][1]==1)) 
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++)
	    udata[(k*x1len + j)*x0len + i] = 0.0;

    if (rank > 1) {
      // x1 left boundary
      if (OnBdry[1][0] && (BdryType[1][0]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x1 right boundary
      if (OnBdry[1][1] && (BdryType[1][1]==1)) 
	for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    if (rank > 2) {
      // x2 left boundary
      if (OnBdry[2][0] && (BdryType[2][0]==1)) 
	for (k=0; k<ugh[2][0]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
      
      // x2 right boundary
      if (OnBdry[2][1] && (BdryType[2][1]==1)) 
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	    for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	      udata[(k*x1len + j)*x0len + i] = 0.0;
    }

  }


//   if (debug)
//     printf("Exiting MFProb::EnforceBoundary routine\n");

  float ftime = MPI_Wtime();
  timers[10] += ftime - stime;

  // return success
  return SUCCESS;

}
#endif
