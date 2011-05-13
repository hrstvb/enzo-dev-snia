/***********************************************************************
/
/  AN EXAMPLE ACTIVE PARTICLE TYPE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EventHooks.h"
#include "ActiveParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_PopIII;

class PopIIIGrid : private grid {
    friend class ActiveParticleType_PopIII;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

class ActiveParticleType_PopIII : public ActiveParticleType
{
public:
  static int EvaluateFormation(grid *thisgrid_orig, float *cooling_time,
			       float *temperature, float *dmfield, 
			       int DensNum, int Vel1Num, int Vel2Num, 
			       int Vel3Num, float DensityUnits);
};

int ActiveParticleType_PopIII::EvaluateFormation
(grid *thisgrid_orig, int DensNum, int Vel1Num, int Vel2Num, int Vel3Num)
{
  SampleParticleGrid *tg =
    static_cast<SampleParticleGrid *>(thisgrid_orig);

  float div;
  int i, j, k, index, offset_y, offset_z;

  float *Refined = tg->BaryonField[tg->NumberOfBaryonFields];

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = tg->GridDimension[0];
  offset_z = tg->GridDimension[0] * tg->GridDimension[1];

  for (k = tg->GridStartIndex[2]; k <= tg->GridEndIndex[2]; k++) {
    for (j = tg->GridStartIndex[1]; j <= tg->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(tg->GridStartIndex[0], j, k);
      for (i = tg->GridStartIndex[0]; i <= tg->GridEndIndex[0]; i++, index++) {

	// 1. Finest level of refinement
	if (Refined[index] != 0.0) continue;
	
	// 2. Density greater than threshold
	if (tg->BaryonField[DensNum] < PopIIIOverDensityThreshold)
	  continue;

	/* 3. Negative divergence: For ZEUS, the velocities are
	   face-centered, and all of the other routines have
	   cell-centered velocities. */

	if (HydroMethod == Zeus_Hydro) {
	  div = tg->BaryonField[Vel1Num][index+1] - tg->BaryonField[Vel1Num][index]
	    + tg->BaryonField[Vel2Num][index+offset_y] - tg->BaryonField[Vel2Num][index]
	    + tg->BaryonField[Vel3Num][index+offset_z] - tg->BaryonField[Vel3Num][index]
	} else {
	  div = tg->BaryonField[Vel1Num][index+1] - 
	    tg->BaryonField[Vel1Num][index-1] + 
	    tg->BaryonField[Vel2Num][index+offset_y] - 
	    tg->BaryonField[Vel2Num][index-offset_y] + 
	    tg->BaryonField[Vel3Num][index+offset_z] - 
	    tg->BaryonField[Vel3Num][index-offset_z]
	}

	if (div > 0.0) continue;

	// 4. t_cool < t_freefall (skip if T > 11000 K)
	dtot = ( tg->BaryonField[DensNum][index] + 
		 dmfield[index] ) * DensityUnits;
	tdyn = sqrt(3.0 * M_PI / 32.0 / GravConst / dtot) / TimeUnits;
	if (tdyn < cooltime[index] && temperature[index] > 1.1e4)
	  continue;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k


  return 0;
}

namespace {
    ActiveParticleType_info *SampleInfo = new ActiveParticleType_info(
            "PopIII", (&ActiveParticleType_PopIII::EvaluateFormation));
}
