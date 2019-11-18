/***********************************************************************
/
/  GRID CLASS (COMPUTE 3D MHD FLUXES)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "../myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int MHDSweepX(float **Prim,  float **Flux3D, int GridDimension[],
	      int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback, FLOAT cellWidth);
int MHDSweepY(float **Prim,  float **Flux3D, int GridDimension[],
	      int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback, FLOAT cellWidth);
int MHDSweepZ(float **Prim,  float **Flux3D, int GridDimension[],
	      int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback, FLOAT cellWidth);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::MHD3D(float **Prim, float **dU, float dt,
		fluxes *SubgridFluxes[], int NumberOfSubgrids,
		float fluxcoef, int fallback)
  /*
     Input:  U[NEQ_SRHYDRO][GridDimension^3]: the conserved variables vector
                                              including ghost zones.
             Prim[NEQ_SRHYDRO+1][GridDimension^3]:
     Output: dU[NEQ_SRHYDRO][activesize^3]: the spatial differenced value without
                                            ghost zones.
             SubgridFluxes:
  */
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
  }

  float *Flux3D[NEQ_MHD+NSpecies+NColor];

  int fluxsize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    fluxsize *= (GridDimension[dim]-2*NumberOfGhostZones+1);
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int Xactivesize = GridDimension[0]-2*NumberOfGhostZones;
  int Yactivesize = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  int Zactivesize = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    Flux3D[field] = new float[fluxsize];
  }
  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    for (int i = 0; i < fluxsize; i++) {
      Flux3D[field][i] = 0.0;
    }
//	  arr_newset(Flux3D + field, fluxsize, 0);
  }

  FLOAT a = 1, dadt;
  FLOAT cellWidth = CellWidth[0][0];
  if(UseBurning && BurningDiffusionMethod==3)
  {
		float *rhoField = Prim[0];
		float *niField = Prim[NEQ_MHD];
		for(size_t i = 0; i < GetGridSize(); i++)
		{
			float rho = rhoField[i];
			if(rho > tiny_number)
			{
				float ni = niField[i];
				niField[i] = (ni > 0) ? (ni / (4 * rho - 3 * ni)) : 0;

			}
			else
			{
				niField[i] = 0;
			}
		}
  }

  /* If using comoving coordinates, multiply dx by a(n+1/2).
     In one fell swoop, this recasts the equations solved by solver
     in comoving form (except for the expansion terms which are taken
     care of elsewhere). */

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
  }

  const int offset[3] = {1, Xactivesize+1, (Xactivesize+1)*(Yactivesize+1)};

  // compute flux at cell faces in x direction
  FLOAT dtdx = dt/(a*CellWidth[0][0]);
  if (MHDSweepX(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback, cellWidth)
      == FAIL) {
    return FAIL;
  }

  // Update dU
  int iflux, ifluxp1, n;
  if (Coordinate == Cartesian) {
    for (int field = 0; field < NEQ_MHD + NSpecies + NColor; field++) {
      n = 0;
      for (int k = 0; k < Zactivesize; k++) {
	for (int j = 0; j < Yactivesize; j++) {
	  iflux = (Xactivesize+1) * (j + k*(Yactivesize+1));
	  for (int i = 0; i < Xactivesize; i++, n++, iflux++) {
	    ifluxp1 = iflux + offset[0];

	    dU[field][n] = -(Flux3D[field][ifluxp1] - Flux3D[field][iflux]) * dtdx;

	    if(field==9 && k == Zactivesize/2 && j==Yactivesize/2)
	    {
//	    	TRACEF("  x-sweep  %lld:  %e  %e", field, Flux3D[field][iflux], Flux3D[field][ifluxp1]);
	    }
	  }
	}
      }
    }

  } // ENDIF Cartesian

  FLOAT geofacr, geofacl;

  if (Coordinate == Cylindrical) {
    FLOAT xl, xr, xc;
    int n = 0, i1;
    for (int k = 0; k < Zactivesize; k++) {
      for (int j = 0; j < Yactivesize; j++) {
	for (int i = 0; i < Xactivesize; i++, n++) {
	  iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
	  i1 = i + GridStartIndex[0];
	  xl = CellLeftEdge[0][i1];
	  xc = xl + 0.5*CellWidth[0][i1];
	  xr = xl + CellWidth[0][i1];
	  geofacr = xr/xc;
	  geofacl = xl/xc;
	  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
            dU[field][n] =
	      -(geofacr * Flux3D[field][iflux+1] - geofacl * Flux3D[field][iflux]) * dtdx;
          }
        }
      }
    }
  }


  if (FluxCorrection) {
    if (this->SaveMHDSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				   Flux3D, 0, fluxcoef, dt) == FAIL) {
      return FAIL;
    }
  }

  // Sweep Y

  if (0 && GridRank > 1) {
    dtdx = dt/(a*CellWidth[1][0]);
    // compute flux in y direction
    if (MHDSweepY(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback, cellWidth)
	== FAIL)
      return FAIL;

    // Update dU
    int n;
    if (Coordinate == Cartesian || Coordinate == Cylindrical) {
      for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
	n = 0;
	for (int k = 0; k < Zactivesize; k++) {
	  for (int j = 0; j < Yactivesize; j++) {
	    iflux = (Xactivesize + 1) * (j+k*(Yactivesize + 1));
	    for (int i = 0; i < Xactivesize; i++, n++, iflux++) {
	      ifluxp1 = iflux + offset[1];
	      dU[field][n] -= (Flux3D[field][ifluxp1]-Flux3D[field][iflux])*dtdx;
//		    if(field==9 && k == Zactivesize/2 && j==Yactivesize/2)
//		    {
//		    	TRACEF("  y-sweep  %lld:  %e  %e", field, Flux3D[field][iflux], Flux3D[field][ifluxp1]);
//		    }
	    }
	  }
	}
      }


      if (FluxCorrection) {
	if (this->SaveMHDSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				       Flux3D, 1, fluxcoef, dt) == FAIL) {
	  return FAIL;
	}
      }
    }
  }  // ENDIF GridRank > 1

  // Sweep Z

  if (0 && GridRank > 2) {
    dtdx = dt/(a*CellWidth[2][0]);
    // compute flux in z direction
    if (MHDSweepZ(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback, cellWidth)
	== FAIL)
      return FAIL;

    // Update dU
    if (Coordinate == Cartesian) {
      int n;
      for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
	n = 0;
	for (int k = 0; k < Zactivesize; k++) {
	  for (int j = 0; j < Yactivesize; j++) {
	    iflux = (Xactivesize+1) * (j + k*(Yactivesize+1));
	    for (int i = 0; i < Xactivesize; i++, n++, iflux++) {
	      ifluxp1 = iflux + offset[2];
	      dU[field][n] -= (Flux3D[field][ifluxp1]-Flux3D[field][iflux])*dtdx;
//    	    if(field==9 && k == Zactivesize/2 && j==Yactivesize/2)
//    	    {
//    	    	TRACEF("  z-sweep  %lld:  %e  %e", field, Flux3D[field][iflux], Flux3D[field][ifluxp1]);
//    	    }
	    }
	  }
	}
      }

    }

    if (Coordinate == Cylindrical) {
      FLOAT xc;
      int n = 0, i1;
      for (int k = 0; k < Zactivesize; k++) {
	for (int j = 0; j < Yactivesize; j++) {
	  for (int i = 0; i < Xactivesize; i++, n++) {
	    iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
	    ifluxp1 = i + (Xactivesize + 1)*(j+(k + 1)*(Yactivesize + 1));
	    i1 = i + GridStartIndex[0];
	    xc = CellLeftEdge[0][i1] + 0.5*CellWidth[0][i1];
	    for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
              dU[field][n] -= (Flux3D[field][ifluxp1] - Flux3D[field][iflux])*dtdx/xc;
            }
          }
        }
      }
    }

    if (FluxCorrection) {
      if (this->SaveMHDSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				     Flux3D, 2, fluxcoef, dt) == FAIL) {
	return FAIL;
      }
    }
  }

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] Flux3D[field];
  }

  return SUCCESS;

}
