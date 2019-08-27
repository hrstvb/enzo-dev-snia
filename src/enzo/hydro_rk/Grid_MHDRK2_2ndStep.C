/***********************************************************************
 /
 /  GRID CLASS (RUNGE-KUTTA SECOND STEP)
 /
 /  written by: Peng Wang
 /  date:       June, 2007
 /  modified1:
 /
 /
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

double ReturnWallTime();
int MHDTimeUpdate_CUDA(float **Prim, int GridDimension[],
int GridStartIndex[], int GridEndIndex[], int GridRank,
float dtdx, float dt, float C_h, float C_p, float cTheta_Limiter);

int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[], float dy[],
float dz[], FLOAT** CellLeftEdges);

int grid::MHDRK2_2ndStep(fluxes *SubgridFluxes[],
int NumberOfSubgrids, int level, ExternalBoundary *Exterior)
/*
 NumberOfSubgrids: the actual number of subgrids + 1
 SubgridFluxes[NumberOfSubgrids]
 */
{

	if(ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}

	if(NumberOfBaryonFields == 0)
	{
		return SUCCESS;
	}TIMER_START("MHDRK2");

#ifdef ECUDA
  if (UseCUDA) {
    this->CudaMHDRK2_2ndStep(SubgridFluxes, NumberOfSubgrids, level, Exterior);
    TIMER_STOP("MHDRK2");
    return SUCCESS;
  }
#endif

	double time1 = ReturnWallTime();

	float *Prim[NEQ_MHD + NSpecies + NColor];
	float *OldPrim[NEQ_MHD + NSpecies + NColor];

	this->ReturnHydroRKPointers(Prim, false);
	this->ReturnOldHydroRKPointers(OldPrim, false);

	if(StellarWindFeedback)
		this->ReduceWindBoundary();

	float *dU[NEQ_MHD + NSpecies + NColor];
	int size = 1;
	for(int dim = 0; dim < GridRank; dim++)
		size *= GridDimension[dim];

	int activesize = 1;
	for(int dim = 0; dim < GridRank; dim++)
		activesize *= (GridDimension[dim] - 2 * NumberOfGhostZones);

	for(int field = 0; field < NEQ_MHD + NSpecies + NColor; field++)
	{
		dU[field] = new float[activesize];
	}

	this->ReturnHydroRKPointers(Prim, true);  //##### added! because Hydro3D needs fractions for species

	/* Compute dU */

	if(OuterVelocitiesClearAtRKStep2Begin)
		ClearOuterVelocities();
//	ClearOuterVelocities(Prim[1], Prim[2], Prim[3], GridDimension[0], GridDimension[1], GridDimension[2], GridRank,
//							CellWidth[0], CellWidth[1], CellWidth[2], CellLeftEdge);

	int fallback = 0;
	if(this->MHD3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 0.5, fallback) == FAIL)
	{
		return FAIL;
	}

	/* Add source terms */

	this->MHDSourceTerms(dU);

	/* Update primitive variables */

	if(this->UpdateMHDPrim(dU, 0.5, 0.5, "MHDRK2_1stStep") == FAIL)
	{
		// fall back to zero order scheme
		fprintf(stderr, "Grid_MHDRK2_2ndStep: Falling back to zero order at RK 2nd step\n");

		this->CopyOldBaryonFieldToBaryonField();
		// change species from density to mass fraction
		for(int ns = NEQ_MHD; ns < NEQ_MHD + NSpecies + NColor; ns++)
		{
			for(int n = 0; n < size; n++)
			{
				Prim[ns][n] /= Prim[iden][n];
			}
		}

		this->ZeroFluxes(SubgridFluxes, NumberOfSubgrids);
		fallback = 1;
		if(this->MHD3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 0.5, fallback) == FAIL)
		{
			return FAIL;
		}
		this->MHDSourceTerms(dU);
		if(this->UpdateMHDPrim(dU, 0.5, 0.5, "MHDRK2_1stStep, fallback") == FAIL)
		{
			fprintf(stderr, "Grid_MHDRK2_2ndStep: Fallback failed, give up...\n");
			return FAIL;
		}
		return FAIL;
	}

	for(int field = 0; field < NEQ_MHD + NSpecies + NColor; field++)
	{
		delete[] dU[field];
	}

	if(OuterVelocitiesClearAtRKStep2End)
		ClearOuterVelocities();
//	ClearOuterVelocities(Prim[1], Prim[2], Prim[3], GridDimension[0], GridDimension[1], GridDimension[2], GridRank,
//							CellWidth[0], CellWidth[1], CellWidth[2], CellLeftEdge);

	TIMER_STOP("MHDRK2");

	/* If we're supposed to be outputting on Density, we need to update
	 the current maximum value of that Density. */

	if(OutputOnDensity == 1)
	{
		int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
		for(int i = 0; i < size; i++)
			CurrentMaximumDensity = max(BaryonField[DensNum][i], CurrentMaximumDensity);
	}

	return SUCCESS;

}
