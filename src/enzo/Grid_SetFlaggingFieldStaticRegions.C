/***********************************************************************
/
/  GRID CLASS (SET THE STATIC REGIONS IN THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

// This routine flags all cells which are adjacent (diagonally counts) to
//  an already flagged Cell.  It also removes all flagged cells which are
//  in the boundary.

#include <stdio.h>
#include <stdlib.h>
#include "myenzoutils.h"
#include "DebugMacros.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::SetFlaggingFieldStaticRegions(int level, int &NumberOfFlaggedCells)
{
  /* declarations */

  int i, j, k, index, region, dim;
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  /* error check */

  if (FlaggingField == NULL) {
    ENZO_FAIL("Flagging Field is undefined.\n");
  }

  for (region = 0; region < MAX_STATIC_REGIONS; region++)
  {
		const FLOAT rInner = StaticRefineShellInnerRadius[region];
		const FLOAT rOuter = StaticRefineShellOuterRadius[region];

		if((rInner < 0) || (rInner > rOuter))
			continue;

//		if(NumberOfFlaggedCells < 0)
//			NumberOfFlaggedCells = 0;

//		if(beforeFlagBufferZones ^ (bool) StaticRefineShellWithBuffer[region])
//			continue;

		if(level > StaticRefineShellLevel[region])
			continue;

		const FLOAT* center = StaticRefineShellCenter[region];
		const FLOAT rInnerSquared = square(rInner);
		const FLOAT rOuterSquared = square(rOuter);
		const FLOAT x0 = center[0];
		FLOAT lxyz[MAX_DIMENSION], rxyz[MAX_DIMENSION];
		long lijk[MAX_DIMENSION], rijk[MAX_DIMENSION];

		for(int dim=0; dim<GridRank;dim++)
		{
			lxyz[dim] = center[dim]-rOuter;
			rxyz[dim] = center[dim]+rOuter;
		}
		if(intersectActive(lxyz, rxyz))
			continue;

		get_ijk_index(lijk, lxyz);
		get_ijk_index(rijk, rxyz);
		if(intersectActive(lijk, rijk))
			continue;

		for(int dim=GridRank; dim<MAX_DIMENSION; dim++)
		{
			lijk[dim] = 0;
			rijk[dim] = 0;
		}

		for(int k = lijk[2]; k <= rijk[2]; k++)
		{
			FLOAT zz = square(CELLCENTER(2, k) - center[2]);
			for(int j = lijk[1]; j < rijk[1]; j++)
			{
				FLOAT yy_zz = square(CELLCENTER(1, j) - center[1]) + zz;
				size_t index = ELT(lijk[0], j, k);
				for(int i = lijk[0]; i <= rijk[0]; i++, index++)
				{
					if(FlaggingField[index])
						continue;
					FLOAT rr = square(CELLCENTER(0, i) - x0) + yy_zz;
					if(rr >= rInnerSquared && rr <= rOuterSquared)
					{
						FlaggingField[index] = 1;
//						NumberOfFlaggedCells++;
					}
				}
			}
		}
  }

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Loop over static regions. */

  for (region = 0; region < MAX_STATIC_REGIONS; region++)
    if (StaticRefineRegionLevel[region] >= level) {

      /* Check if there is any overlap. */

      int Overlap = TRUE;
      for (dim = 0; dim < GridRank; dim++) {
	Left[dim] = max(StaticRefineRegionLeftEdge[region][dim],
			GridLeftEdge[dim]);
	Right[dim] = min(StaticRefineRegionRightEdge[region][dim],
			GridRightEdge[dim]);
	if (Left[dim] >= Right[dim])
	  Overlap = FALSE;
      }

      if (Overlap == TRUE) {
	int Start[] = {0,0,0}, End[] = {0,0,0};
	for (dim = 0; dim < GridRank; dim++) {
	  Start[dim] = nint((Left[dim] - CellLeftEdge[dim][0])/
			    CellWidth[dim][0]);
	  End[dim] = nint((Right[dim] - CellLeftEdge[dim][0])/
			  CellWidth[dim][0]) - 1;
	}

	for (k = Start[2]; k <= End[2]; k++)
	  for (j = Start[1]; j <= End[1]; j++) {
	    index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
	    for (i = Start[0]; i <= End[0]; i++, index++)
	      FlaggingField[index] = 1;
	  }

      } // end: if (Overlap)

    } // end: if (StaticRefineRegionLevel[dim] == level)

  /* Count up the number of flagged cells & report. */
  NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    NumberOfFlaggedCells += FlaggingField[i];

  return SUCCESS;
}
