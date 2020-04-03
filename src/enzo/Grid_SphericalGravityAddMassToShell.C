/*
 * SphericalGravityAddMassToShell
 * dcollins.  Oct. 16 2018. 15:15.
 */
#include <stdio.h>
#include "hdf5.h"

#include "myenzoutils.h"
#include "DebugMacros.h"
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "fortran.def"

/* function prototypes */

//int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
//int FindField(int f, int farray[], int n);
int SphericalGravityAllocateBins(size_t** countBins, FLOAT** massBins, FLOAT** binAccels, FLOAT** binAccelSlopes,
FLOAT** centersOfMassBins, FLOAT** kineticEBins, FLOAT** magneticEBins, FLOAT** volumeBins, int domainRank);

int SphericalGravityComputeBinIndex(FLOAT r);

//bool PointInChildren(const FLOAT* const point, const HierarchyEntry* const firstChild)
//{
//	const HierarchyEntry* child = firstChild;
//	while(child)
//	{
//		grid* g = child->GridData;
//		if(g && g->PointInGrid(point))
//			return true;
//		child = child->NextGridThisLevel;
//	}
//	return false;
//}
//
//bool PointInChildren(const FLOAT* const point, const LevelHierarchyEntry* const myLevelHierarchyEntry)
//{
//	const HierarchyEntry* const he = (myLevelHierarchyEntry) ? myLevelHierarchyEntry->GridHierarchyEntry : NULL;
//	return PointInChildren(point, (he) ? he->NextGridNextLevel : NULL);
//}

int grid::SphericalGravityAddMassToShell(size_t* countBins, FLOAT* densBins, FLOAT** cmBins, FLOAT** kinEBins,
FLOAT** magEBins, LevelHierarchyEntry* myLevelHierarchyEntry)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	HierarchyEntry* he = myLevelHierarchyEntry->GridHierarchyEntry->ParentGrid;
	grid* mypg = (he) ? he->GridData : NULL;
	int myLevel = 0;
	for(; he; he = he->ParentGrid)
		myLevel++;
	int maxRefLevel = min(max(0, MaximumRefinementLevel), MAX_DEPTH_OF_HIERARCHY);
	int maxGravityLevel = (SphericalGravityMaxHierarchyLevel < 0) ? maxRefLevel : SphericalGravityMaxHierarchyLevel;
	maxGravityLevel = min(maxGravityLevel, maxRefLevel);
	if(myLevel > maxGravityLevel)
		return SUCCESS;
	bool checkChildren = myLevel < maxGravityLevel;
	const size_t &N = SphericalGravityActualNumberOfBins;

	float *densField, *vFields[3], *BFields[3];

	MHD_SNIA_GetFields(&densField, NULL, NULL, NULL, NULL, NULL, vFields, NULL, NULL, NULL, BFields, NULL, NULL, NULL);

	// Add quantities to the bins for this grid.
	// Cell conunt in each bin is literal.
	// For the mass and the energies components we add up
	// the corresponding densities and multiply by the cell
	// volume after the loop.
	//Loop variables

	HierarchyEntry* firstChild = (myLevelHierarchyEntry) ? myLevelHierarchyEntry->GridHierarchyEntry : NULL;
	firstChild = (firstChild) ? firstChild->NextGridNextLevel : NULL;

	FLOAT dens, r, xVec[3], rVec[3], xx, yy_zz, zz;
	xx = yy_zz = zz = 0;
	arr_set(xVec, 3, 0);
	arr_set(rVec, 3, 0);
	size_t rbin, index;
	switch(GridRank)
	{
	case 1:
		for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
		{
			r = rVec[0] = (xVec[0] = CELLCENTER(0, i)) - SphericalGravityCenter[0];
			// If this point is in any of the children, skip it.
			if(checkChildren)
				if(PointInChildrenActiveNB(xVec, firstChild))
					continue;

			if(-1 == (rbin = SphericalGravityComputeBinIndex(r)))
				continue;
			countBins[rbin]++;
			densBins[rbin] += densField[i]; // * cellVolume;

			FLOAT **bins = cmBins;
			for(int dim = 0; dim < GridRank; bins++, dim++)
				*bins[rbin] += dens * rVec[dim];
			bins = kinEBins;
			for(int dim = 0; dim < GridRank; bins++, dim++)
				*bins[rbin] += dens * square(vFields[dim][index]);
			bins = magEBins;
			for(int dim = 0; dim < 3; bins++, dim++)
				if(*bins)
					*bins[rbin] += square(BFields[dim][index]) / dens;
			index++;
		}
		break;
	case 2:
		for(int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
		{
			yy_zz = square(rVec[1] = (xVec[1] = CELLCENTER(1, j)) - SphericalGravityCenter[1]);
			index = GridDimension[0] * j;
			for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
			{
				xx = square(rVec[0] = (xVec[0] = CELLCENTER(0, i)) - SphericalGravityCenter[0]);
				r = sqrt(xx + yy_zz);
				// If this point is in any of the children, skip it.
				if(checkChildren)
					if(PointInChildrenActiveNB(xVec, firstChild))
						continue;

				if(-1 == (rbin = SphericalGravityComputeBinIndex(r)))
					continue;
				countBins[rbin]++;
				densBins[rbin] += densField[index]; // * cellVolume;

				FLOAT **bins = cmBins;
				for(int dim = 0; dim < GridRank; bins++, dim++)
					*bins[rbin] += dens * rVec[dim];
				bins = kinEBins;
				for(int dim = 0; dim < GridRank; bins++, dim++)
					*bins[rbin] += dens * square(vFields[dim][index]);
				bins = magEBins;
				for(int dim = 0; dim < 3; bins++, dim++)
					if(*bins)
						*bins[rbin] += square(BFields[dim][index]) / dens;
				index++;
			}
		}
		break;
	default:
		bool printdebuginfo = true;
		for(int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
		{
			zz = square(rVec[2] = (xVec[2] = CELLCENTER(2, k)) - SphericalGravityCenter[2]);
			for(int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
			{
				yy_zz = square(rVec[1] = (xVec[1] = CELLCENTER(1, j)) - SphericalGravityCenter[1]) + zz;
				index = GridDimension[0] * (j + GridDimension[1] * k);
				for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
				{
					xx = square(rVec[0] = (xVec[0] = CELLCENTER(0, i)) - SphericalGravityCenter[0]);
					r = sqrt(xx + yy_zz);

					// If this point is in any of the children, skip it.
					if(checkChildren)
						if(PointInChildrenActiveNB(xVec, firstChild))
							continue;

//					if(0 && j == 55 && k == 55)
//						TRACEGF(" mylevel=%lld   %lld %lld %lld   %e  %e  %e", myLevel, i, j, k, xVec[0], xVec[1],
//								xVec[2]);

					if(-1 == (rbin = SphericalGravityComputeBinIndex(r)))
						continue;

//					if(j == (GridDimension[1] / 2) && k == (GridDimension[2] / 2) && i <= GridDimension[0] / 2)
//						printf("i,j,k=%03d,%04d,%04d, (%4f,%4f,%4f), r=%4f, rho=%e, bin=%lld\n", i, j, k,
//								rVec[0] * 1e-5, rVec[1] * 1e-5, rVec[2] * 1e-5, r * 1e-5, dens * 1e-9, rbin);

					densBins[rbin] += (dens = densField[index]); // * cellVolume;

					//The other bins are optional
					if(countBins)
						countBins[rbin]++;
					if(cmBins && cmBins[0])
						for(int dim = 0; dim < GridRank; dim++)
							cmBins[dim][rbin] += dens * rVec[dim];
					if(kinEBins && kinEBins[0])
						for(int dim = 0; dim < GridRank; dim++)
							kinEBins[dim][rbin] += dens * square(vFields[dim][index]);
					if(magEBins && magEBins[0])
						for(int dim = 0; dim < GridRank; dim++)
							magEBins[dim][rbin] += square(BFields[dim][index]) / dens;

					index++;
				} // for i, x, dim0
			} // for j, y, dim1
		} // for k, z, dim2
	} //switch(GridRank)

//Add the grid bins to the bins that are global on this processor.
	FLOAT cellVolume = getCellVolume();
//	{ //debug
//		size_t count = arr_sum(countBins, N) / ((mypg) ? 8 : 1);
//		int mypgid = (mypg) ? mypg->GetGridID() : -1;
//		arr_ax(countBins, N, (mypg == NULL) ? 8 : 1);
//		TRACEGF("%lld  %lld  %e  %lld %lld", count, GetActiveSize(), cellVolume, ID, mypgid);
//	}

//	if(MyProcessorNumber == ROOT_PROCESSOR)
//	{
//		char* sbuf = new char[64 * N];
//		char* s = sbuf;
//		s += sprintf(s, "SphericalGravityActualNumberOfBins = %lld\n", N);
//		s += sprintf(s, "cellVolume = %e\n", cellVolume);
//		s += sprintfvec(s, "SphericalGravityShellMasses = {", "%lld:%e", ", ", "}\n", SphericalGravityShellMasses, N,
//						true, true);
//		s += sprintfvec(s, "SphericalGravityInteriorMasses = {", "%lld:%e", ", ", "}\n", densBins,
//						N, true, true);
//		printf(sbuf);
//		delete sbuf;
//	}

//	if(SphericalGravityShellCellCounts && countBins)
//		arr_xpy(SphericalGravityShellCellCounts, countBins, N);
//	if(SphericalGravityShellVolumes && countBins)
//		arr_axpy(SphericalGravityShellVolumes, countBins, N, cellVolume);
	arr_axpy(SphericalGravityShellMasses, densBins, N, cellVolume);
//	for(int dim = 0; dim < GridRank; dim++)
//		arr_axpy(SphericalGravityShellCentersOfMass[dim], cmBins[dim], N, cellVolume);
//	for(int dim = 0; dim < GridRank; dim++)
//		arr_axpy(SphericalGravityShellKineticEnergies[dim], kinEBins[dim], N, 0.5 * cellVolume);
//	for(int dim = 0; dim < 3; dim++)
//		arr_axpy(SphericalGravityShellMagneticEnergies[dim], magEBins[dim], N, 0.5 * cellVolume);

	return SUCCESS;
}

int grid::SphericalGravityAddMassToShell()
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

// Allocate bins for this grid.
	size_t* countBins = NULL;
	FLOAT* densBins = NULL;
	FLOAT** cmBins = arr_newset<FLOAT*>(MAX_DIMENSION, NULL);
	FLOAT** kinEBins = arr_newset<FLOAT*>(MAX_DIMENSION, NULL);
	FLOAT** magEBins = (UseMHD || UseMHDCT) ? arr_newset<FLOAT*>(MAX_DIMENSION, NULL) : NULL;

	SphericalGravityAllocateBins(&countBins, &densBins, NULL, NULL, cmBins, kinEBins, magEBins, NULL, GridRank);
	int retval = SphericalGravityAddMassToShell(countBins, densBins, cmBins, kinEBins, magEBins, NULL);

	delete countBins;
	delete densBins;
	delete[] cmBins;
	delete[] kinEBins;
	delete[] magEBins;
	return retval;
	return SUCCESS;
}
