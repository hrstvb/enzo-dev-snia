/***********************************************************************
/
/  GRID CLASS (IS THIS PHOTON CROSSING AN OPTICALLY THIN CELL)
/
/  written by: John Wise
/  date:       December, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

/* If another ray from the same source has already touched this cell,
   skip computation and go to the next cell.  Estimate flux at the
   cell center, while we're at it. */

int grid::OpticallyThinPhoton(PhotonPackageEntry **PP, float min_tau, int index,
			      FLOAT radius, float dx, float dx2, FLOAT ce[],
			      FLOAT &scaling, bool &skip, bool &thin)
{

  int dim;
  FLOAT cr2, dxc;

  if ((*PP)->ColumnDensity < min_tau) {

    thin = true;

    if ((RayMarker[index] >> (*PP)->SourceNumber & 1) == 0) {
      // Radius from source to cell center
      for (dim = 0, cr2 = 0.0; dim < MAX_DIMENSION; dim++) {
	dxc = (*PP)->SourcePosition[dim] - (ce[dim] + 0.5 * dx);
	cr2 += dxc * dxc;
      }
      cr2 = max(cr2, 0.01*dx2); // make sure the rates don't go infinite
      scaling = radius*radius / cr2;
      skip = false;
    } // ENDIF RayMarker == false
    else
      // RayMarker == true -> propagate ray but no contribute to rates
      skip = true;

  } // ENDIF OpticallyThin
  else {
    if ((RayMarker[index] >> (*PP)->SourceNumber & 1) == 1) {
      skip = true;
      thin = true;
    } else {
      skip = false;
      thin = false;
    }
  }

  return SUCCESS;

}

/* Geometric correction factor because the ray's solid angle could not
   completely cover the cell */

float RayGeometricCorrection(const FLOAT oldr, 
			     const FLOAT radius,
			     const FLOAT ddr, 
			     const FLOAT s[], 
			     const float u[],
			     const FLOAT ce[],
			     const float dxhalf,
			     const float dtheta,
			     const FLOAT roundoff)
{

  int dim;
  float midpoint, nearest_edge;
  float m[3], slice_factor, slice_factor2, sangle_inv;

  midpoint = oldr + 0.5f*ddr - roundoff;
  nearest_edge = -1e20;
  for (dim = 0; dim < 3; dim++)
    m[dim] = fabs(s[dim] + midpoint * u[dim] - (ce[dim] + dxhalf));
  nearest_edge = max(max(m[0], m[1]), m[2]);
  sangle_inv = 1.0 / (dtheta*radius);
  slice_factor = min(0.5f + (dxhalf-nearest_edge) * sangle_inv, 1.0f);
  slice_factor2 = slice_factor * slice_factor;
  //slice_factor2 = 1.0f;

  return slice_factor2;

}
