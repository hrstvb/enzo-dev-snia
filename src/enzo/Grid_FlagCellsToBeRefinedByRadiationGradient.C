/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE OF RADIATION FIELD)
/
/  written by: Daniel R. Reynolds
/  date:       September, 2011
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::FlagCellsToBeRefinedByRadiationGradient()
{
  // declarations 
  int i, j, k, index, dim;
 
  // Return if this grid is not on this processor
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  // error check
  if (FlaggingField == NULL)
    ENZO_FAIL("Flagging Field is undefined.");
 
  // Make sure quantities are defined at least to dim 3 (if GridRank < 3)
  for (dim=GridRank; dim<3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  // compute total grid size 
  int size=1;
  for (dim=0; dim<GridRank; dim++)   size *= GridDimension[dim];
 
  // initialize return value
  int NumberOfFlaggedCells = 0;
 
  // only take gradient over active dimensions
  int xOffset = (GridDimension[0] > 1) ? 1 : 0;
  int yOffset = (GridDimension[1] > 1) ? GridDimension[0]+1 : 0;
  int zOffset = (GridDimension[2] > 1) ? GridDimension[0]*GridDimension[1]+1 : 0;
 
  // get Radiation field index 
  float *RadField = this->AccessRadiationFrequency0();
  if (RadField == NULL)
    ENZO_FAIL("Could not find RadiationFrequency0 field.");

  // use first entry in MinimumSlopeForRefinement array to determine refinement 
  // threshold; if left unset, just return
  float MinRefinementSlope = MinimumSlopeForRefinement[0];
  if (MinRefinementSlope == 0.0)  return NumberOfFlaggedCells;

  // iterate over the grid, marking cells as needed
  float maxgrad=0.0, mingrad=0.0, avggrad=0.0;
  float atol = 0.1;
  float gradx, grady, gradz, gradnorm;
  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
	index = i + (j + k*GridDimension[1])*GridDimension[0];
	gradx = (RadField[index+xOffset] - RadField[index-xOffset]) /
	  max(fabs(RadField[index]), atol);
	grady = (RadField[index+yOffset] - RadField[index-yOffset]) /
	  max(fabs(RadField[index]), atol);
	gradz = (RadField[index+zOffset] - RadField[index-zOffset]) /
	  max(fabs(RadField[index]), atol);
	gradnorm = sqrt(gradx*gradx + grady*grady + gradz*gradz);
	maxgrad = (gradnorm > maxgrad) ? gradnorm : maxgrad;
	mingrad = (gradnorm < mingrad) ? gradnorm : maxgrad;
	avggrad += gradnorm;
	if (gradnorm > MinRefinementSlope) {
	  FlaggingField[index] += 1;
	  NumberOfFlaggedCells++;
	}
      }  // end loop over cells
  
  avggrad /= ((GridEndIndex[2]-GridStartIndex[2]+1) *
	      (GridEndIndex[1]-GridStartIndex[1]+1) * 
	      (GridEndIndex[0]-GridStartIndex[0]+1) );
  printf("FlagCellsRadiation: grad max = %g, min = %g, avg = %g, NFlagCells = %i\n",maxgrad,mingrad,avggrad,NumberOfFlaggedCells);
	  
  return NumberOfFlaggedCells;
}
