/***********************************************************************
/
/  GRID CLASS (DEFRAGMENTS BAYON MEMORY POOL)
/
/  written by: Tom Abel
/  date:       September 2011
/  modified1:
/
/  PURPOSE:
/       Copy all Baryon fields to a new baryon pool getting rid of wasted 
/       space.
/       Should only be called from DefragmentMemoryPools.C
/ 
************************************************************************/
 
//  Compute derived quantites
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::DefragmentBaryonMemoryPool()
{
 
  /*  Baryons only: set up field quantities and allocate fields
      (we assume here the grid is uniform in each dimension) */
 

  FLOAT AllCellWidth, GridLeftIncludingBoundary;
  int size=1, gsize=1;
  
  for (int dim=0; dim<GridRank; dim++) 
    {
      size *= GridDimension[dim];
      gsize *= GravitatingMassFieldDimension[dim];
    }

  for (int field = 0; field <NumberOfBaryonFields; field++) 
    {
      if (BaryonField[field] != NULL)
	{
	  float *temp = AllocateNewBaryonField(size);  // new space 
	  for (int index=0; index < size; index++)
	    temp[index] = BaryonField[field][index];   // copy correct values
	  FreeBaryonFieldMemory(BaryonField[field]);
	  BaryonField[field] = temp;                   // point BaryonField to correct position
	} //  note we did not free the old space. We do that in DefragmentMemoryPools.C
    }
 
//  OldBaryonFields 
  for (int field = 0; field < NumberOfBaryonFields; field++) 
    {
      if (OldBaryonField[field] != NULL) 
	{
	  float *temp = AllocateNewBaryonField(size);  // new space 
	  for (int index=0; index < size; index++)
	    temp[index] = OldBaryonField[field][index];   // copy correct values
	  FreeBaryonFieldMemory(OldBaryonField[field]);
	  OldBaryonField[field] = temp;    
	}
    }
 
  // if (BaryonField[NumberOfBaryonFields] != NULL) 
  //   {
  //   float *temp = AllocateNewBaryonField(size);  // new space 
  //   for (int index=0; index < size; index++) 
  //     temp[index] = BaryonField[NumberOfBaryonFields][index];   // copy correct values
  //   FreeBaryonFieldMemory(BaryonField[NumberOfBaryonFields]);
  //   BaryonField[NumberOfBaryonFields] = temp;                   // point BaryonField to correct position
  // }

  if (GravitatingMassField != NULL) 
    {
    float *temp = AllocateNewBaryonField(gsize);  // new space 
    for (int index=0; index < gsize; index++)
      temp[index] = GravitatingMassField[index];   // copy correct values
    FreeBaryonFieldMemory(GravitatingMassField);
    GravitatingMassField = temp;                   // point BaryonField to correct position
    }

  if (GravitatingMassFieldParticles != NULL) 
    {
    float *temp = AllocateNewBaryonField(gsize);  // new space 
    for (int index=0; index < gsize; index++)
      temp[index] = GravitatingMassFieldParticles[index];   // copy correct values
    FreeBaryonFieldMemory(GravitatingMassFieldParticles);
    GravitatingMassFieldParticles = temp;                   
    }

  if (PotentialField != NULL) 
    {
    float *temp = AllocateNewBaryonField(gsize);  // new space 
    for (int index=0; index < gsize; index++)
      temp[index] = PotentialField[index];   // copy correct values
    FreeBaryonFieldMemory(PotentialField);
    PotentialField = temp;
    }


  //  FreeBaryonFieldMemory(AccelerationField[dim]);
// #ifdef SAB
//     FreeBaryonFieldMemory(OldAccelerationField[dim]) = NULL;
// #endif

  if (BoundaryFluxes != NULL)
    {
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	int bfsize = 1;
	for (int i = 0; i < GridRank; i++)
	  bfsize *= BoundaryFluxes->LeftFluxEndGlobalIndex[dim][i] -
	    BoundaryFluxes->LeftFluxStartGlobalIndex[dim][i] + 1;
	for (int field = 0; field < NumberOfBaryonFields; field++) {
	  float *templ = AllocateNewBaryonField(bfsize);
	  float *tempr = AllocateNewBaryonField(bfsize);
	  for (int index=0; index< bfsize; index++) 
	    {
	      templ[index] = BoundaryFluxes->LeftFluxes[field][dim][index];
	      tempr[index] = BoundaryFluxes->RightFluxes[field][dim][index];
	    };
	  FreeBaryonFieldMemory(BoundaryFluxes->LeftFluxes[field][dim]);
	  FreeBaryonFieldMemory(BoundaryFluxes->RightFluxes[field][dim]);
    	  BoundaryFluxes->LeftFluxes[field][dim]  = templ;
	  BoundaryFluxes->RightFluxes[field][dim] = tempr;
	}
      } // end: loop over dims
    }
  return;
}





 

 
