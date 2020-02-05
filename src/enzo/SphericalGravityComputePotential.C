/***********
 * ComputeSphericalGravityPotential
 * dcollins.  October 16 2018.  14:58.
 * Radially bins mass into shells and interior mass.
 * *********/
#include "preincludes.h"

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "myenzoutils.h"
#include "MHDInitialProfile.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "LevelArrayIterator.h"

#include "DebugMacros.h"

int CommunicationBroadcastValues(FLOAT *Values, int Number, int BroadcastProcessor);

int SphericalGravityAllocateBins(size_t** countBins, FLOAT** massBins, FLOAT** binAccels,
FLOAT** binAccelSlopes, FLOAT** centersOfMassBins, FLOAT** kineticEBins, FLOAT** magneticEBins, FLOAT** volumeBins,
int domainRank)
{
	size_t &N = SphericalGravityActualNumberOfBins;
	if(countBins)
		arr_delnewset(countBins, N, 0);
	if(massBins)
		arr_delbrnewset(massBins, N, 0);
	if(binAccels)
		arr_delbrnewset(binAccels, N, 0);
	if(binAccelSlopes)
		arr_delbrnewset(binAccelSlopes, N, 0);
	if(centersOfMassBins)
		for(int dim = 0; dim < domainRank; dim++)
			arr_delnewset(centersOfMassBins + dim, N, 0);
	if(kineticEBins)
		for(int dim = 0; dim < domainRank; dim++)
			arr_delbrnewset(kineticEBins + dim, N, 0);
	if(magneticEBins && (UseMHD || UseMHDCT))
		for(int dim = 0; dim < domainRank; dim++)
			arr_delbrnewset(magneticEBins + dim, N, 0);
	if(volumeBins)
		arr_delnewset(volumeBins, N, 0);

	return SUCCESS;
}

int SphericalGravityClearBins(size_t* countBins, FLOAT* massBins, FLOAT** binAccels,
FLOAT** binAccelSlopes, FLOAT** centersOfMassBins, FLOAT** kineticEBins, FLOAT** magneticEBins, FLOAT** volumeBins,
int domainRank)
{
	size_t &N = SphericalGravityActualNumberOfBins;
	if(countBins)
		arr_set(countBins, N, 0);
	if(massBins)
		arr_set(massBins, N, 0);
	if(binAccels)
		arr_set(binAccels, N, 0);
	if(binAccelSlopes)
		arr_set(binAccelSlopes, N, 0);
	if(centersOfMassBins)
		for(int dim = 0; dim < domainRank; dim++)
			arr_set(centersOfMassBins[dim], N, 0);
	if(kineticEBins)
		for(int dim = 0; dim < domainRank; dim++)
			arr_set(kineticEBins[dim], N, 0);
	if(magneticEBins && (UseMHD || UseMHDCT))
		for(int dim = 0; dim < domainRank; dim++)
			arr_set(magneticEBins[dim], N, 0);
	if(volumeBins)
		arr_set(volumeBins, N, 0);

	return SUCCESS;
}

int SphericalGravityDetermineUniformBins()
{
//	if(UseSpherGrav)
//		return SpherGravDetermineUniformBins();

	int nOk = SphericalGravityInnerRadius >= 0;
	nOk += SphericalGravityNumberOfBins > 0;
	nOk += SphericalGravityBinSize > 0;
	nOk += SphericalGravityOuterRadius > max(0, SphericalGravityInnerRadius);
	if(nOk < 3)
	{
		fprintf(stderr, "FAILURE: At least 3 of the following conditions need to be met for uniform bins:"
				"SphericalGravityBinSize > 0, SphericalGravityNumberOfBins > 0, "
				"SphericalGravityInnerRadius >= 0, SphericalGravityOuterRadius > max(0, SphericalGravityInnerRadius),"
				"but the values are %"FSYM", %"ISYM", %"FSYM", %"FSYM".\n",
				SphericalGravityBinSize, SphericalGravityNumberOfBins, SphericalGravityInnerRadius,
				SphericalGravityOuterRadius);
		ENZO_FAIL("SphericalGravityInitializeUniformBins: "
					"Improper parameters inner/outer radius, bin number and/or bin size.\n");
	}

	if(nOk == 4)
	{
//		printf("SphericalGravityInitializeUniformBins: Warning: "
//				"Ignored parameter SphericalGravityOuterRadius. "
//				"Uniform bis caalculated based on the inner radius, "
//				"the bin number and size.\n");
	}

	if(nOk == 4 || SphericalGravityOuterRadius < 0)
	{
		//nothing to do.
	}
	else if(SphericalGravityInnerRadius < 0)
	{
		SphericalGravityInnerRadius = SphericalGravityOuterRadius
				- SphericalGravityNumberOfBins * SphericalGravityBinSize;
		if(SphericalGravityInnerRadius < 0)
			ENZO_FAIL("SphericalGravityInitializeUniformBins: "
						"Improper values fo outer radius, bin number, bin size."
						"The calculated inner radius is negative.\n");
	}
	else
	{
		if(SphericalGravityInnerRadius >= SphericalGravityOuterRadius)
			ENZO_FAIL("SphericalGravityInitializeUniformBins: "
						"Improper parameters inner and outer radii."
						"The inner radius should be < the outer radius.\n");

		if(SphericalGravityBinSize <= 0)
		{
			SphericalGravityBinSize = (SphericalGravityOuterRadius - SphericalGravityInnerRadius)
					/ SphericalGravityNumberOfBins;
		}
		if(SphericalGravityNumberOfBins <= 0)
		{
			SphericalGravityNumberOfBins = int(
					(SphericalGravityOuterRadius - SphericalGravityInnerRadius) / SphericalGravityBinSize);
			SphericalGravityBinSize = (SphericalGravityOuterRadius - SphericalGravityInnerRadius)
					/ SphericalGravityNumberOfBins;
		}
	}

	SphericalGravityHasCentralBin = SphericalGravityInnerRadius > 0;
	size_t &N = SphericalGravityActualNumberOfBins;
	N = SphericalGravityNumberOfBins + SphericalGravityHasCentralBin + 1;
//	TRACEF("%lld x %e = %e .. %e", N, SphericalGravityActualNumberOfBins, SphericalGravityBinSize,
//			SphericalGravityInnerRadius, SphericalGravityOuterRadius);

	arr_delnewset(&SphericalGravityBinLeftEdges, N, SphericalGravityBinSize);

	SphericalGravityBinLeftEdges[0] = 0.0;
	if(SphericalGravityHasCentralBin)
		SphericalGravityBinLeftEdges[1] = SphericalGravityInnerRadius;
	arr_cumsum(SphericalGravityBinLeftEdges, N, 0);
//	for (int i = 1; i < SphericalGravityActualNumberOfBins - 1; i++)
//		SphericalGravityBinRightEdges[i] = SphericalGravityBinLeftEdges[i + 1];
//	for (int i = 0; i < SphericalGravityActualNumberOfBins - 1; i++)
//		SphericalGravityBinCenters[i] = 0.5 * (SphericalGravityBinLeftEdges[i] + SphericalGravityBinRightEdges[i]);

	return SUCCESS;
}

int SphericalGravityDetermineBins()
{
	if(SphericalGravityUniformBins)
	{
		if(FAIL == SphericalGravityDetermineUniformBins())
			ENZO_FAIL("Non-uniform bins for spherical gravity are not implemented.\n");
	}
	else
	{		//	Initialize non-uniform bins here. This may require communication.
		ENZO_FAIL("Non-uniform bins for spherical gravity are not implemented.\n");
	}

	return SUCCESS;
}

size_t SphericalGravityComputeUniformBinIndex(FLOAT r)
{
	if(r < SphericalGravityInnerRadius)
		return 0;
	if(r > SphericalGravityOuterRadius)
		return SphericalGravityActualNumberOfBins - 1;
	size_t i = (r - SphericalGravityInnerRadius) / SphericalGravityBinSize + SphericalGravityHasCentralBin;
	return i;
}

size_t SphericalGravityComputeBinIndex(FLOAT r)
{
	if(SphericalGravityUniformBins)
		return SphericalGravityComputeUniformBinIndex(r);
	size_t i = findmaxlte(SphericalGravityBinLeftEdges, SphericalGravityActualNumberOfBins, r);
	return i;
}

//int SpherGravComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData)
//{
//	if(SpherGravDetermineBins() == FAIL)
//		ENZO_FAIL("Coudn't determine spherical gravity bins.\n");
//
//	int L = SphericalGravityActualNumberOfBins + 1;
//	for(int level = 0; level < L; level++)
//	{
//		SpherGravAllocateBins(level, SpherGravShellCellCounts + level, SpherGravShellMasses + level,
//								SpherGravBinAccels + level, SpherGravBinAccelSlopes + level,
//								SpherGravShellCentersOfMass, SpherGravShellKineticEnergies,
//								SpherGravShellMagneticEnergies, SpherGravShellVolumes + level, MetaData->TopGridRank);
//	}
//
//	size_t maxSubgridCount = getMaxSubgridCount(LevelArray);
//	HierarchyEntry* subgrids = new HierarchyEntry[maxSubgridCount];
//
//	LevelHierarchyEntry* lhe = LevelArray[0];
//	while(lhe)
//	{
//		grid* g = lhe->GridHierarchyEntry->GridData;
//		HierarchyEntry* subg = lhe->GridHierarchyEntry->NextGridNextLevel;
////		int nsubgrids = getSubgrids()subgrids, lhe);
//		HierarchyEntry* subgrids[1024];
//		while(subg)
//		{
////			subgrids[nsubgrids++] = subg;
//			subg = subg->NextGridThisLevel;
//		}
//
//		// get subgrid edges in int coords.
//		g->GetGridLeftEdge();
//		lhe = lhe->NextGridThisLevel;
//	}
//
//	return SUCCESS;
//}

int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData,
	bool interpolateMissingShells)
{
	if(UseSphericalGravity == 0)
		return SUCCESS;

	SphericalGravityMaxHierarchyLevel = min(SphericalGravityMaxHierarchyLevel, MaximumRefinementLevel);
//	UseSpherGrav = SphericalGravityMaxHierarchyLevel > 0;
//	if(UseSpherGrav)
//		return SpherGravCopmutePotential(LevelArray, MetaData);

	if(SphericalGravityDetermineBins() == FAIL)
		ENZO_FAIL("Coudn't determine spherical gravity bins.\n");
	SphericalGravityAllocateBins(&SphericalGravityShellCellCounts, &SphericalGravityShellMasses,
									&SphericalGravityBinAccels, &SphericalGravityBinAccelSlopes,
									SphericalGravityShellCentersOfMass, SphericalGravityShellKineticEnergies,
									SphericalGravityShellMagneticEnergies, &SphericalGravityShellVolumes,
									MetaData->TopGridRank);

	const size_t &N = SphericalGravityActualNumberOfBins;

// Compute potential based on the top level (0) or all levels (1)?
#if(1)
	LevelArrayIterator it = LevelArrayIterator(LevelArray);

	size_t* countBins = NULL;
	FLOAT* densBins = NULL;
	FLOAT** cmBins = arr_newset<FLOAT*>(MAX_DIMENSION, NULL);
	FLOAT** kinEBins = arr_newset<FLOAT*>(MAX_DIMENSION, NULL);
	FLOAT** magEBins = (UseMHD || UseMHDCT) ? arr_newset<FLOAT*>(MAX_DIMENSION, NULL) : NULL;

	SphericalGravityAllocateBins(&countBins, &densBins, NULL, NULL, cmBins, kinEBins, magEBins, NULL, 3);
	int maxRefLevel = max(0, min(MaximumRefinementLevel, MAX_DEPTH_OF_HIERARCHY));
	int maxGravityLevel = min(SphericalGravityMaxHierarchyLevel, maxRefLevel);
	for(grid* g = it.firstFromTop(); g; g = it.next())
	{
		if(it.currentLevel > maxGravityLevel)
			break;
		SphericalGravityClearBins(countBins, densBins, NULL, NULL, cmBins, kinEBins, magEBins, NULL, 3);
		int retval = g->SphericalGravityAddMassToShell(countBins, densBins, cmBins, kinEBins, magEBins,
														it.currentEntry);
	}

	delete countBins;
	delete densBins;
	delete[] cmBins;
	delete[] kinEBins;
	delete[] magEBins;

//	TRACEF("Total count per MPI process: #%lld  %lld ", MyProcessorNumber, arr_sum(SphericalGravityShellCellCounts, N));
//	CommunicationSumValues(SphericalGravityShellCellCounts, N);
//	if(MyProcessorNumber==ROOT_PROCESSOR)
//		TRACEF("Total count for the domain: #%lld  %lld ", MyProcessorNumber, arr_sum(SphericalGravityShellCellCounts, N));

#else
	//	if(Top level gravity only
	LevelHierarchyEntry* lhe = LevelArray[0];
	while(lhe)
	{
		lhe->GridData->SphericalGravityAddMassToShell();
		lhe = lhe->NextGridThisLevel;
	}
#endif

//It might be better to use MPI_Alltoallv; see CommunicationShareParticles.C
//That takes some more setup, so in the words of Mike, "make it work then make it work fast."
	CommunicationSumValues(SphericalGravityShellMasses, N);
//	for (int dim = 0; dim < MAX_DIMENSION; dim++)
//		if (SphericalGravityShellCentersOfMass[dim])
//		{
//			CommunicationSumValues(SphericalGravityShellCentersOfMass[dim], N);
//		}
//	for (int dim = 0; dim < MAX_DIMENSION; dim++)
//		if (SphericalGravityShellKineticEnergies[dim])
//		{
//			CommunicationSumValues(SphericalGravityShellKineticEnergies[dim], N);
//		}
//	for (int dim = 0; dim < MAX_DIMENSION; dim++)
//		if (SphericalGravityShellMagneticEnergies[dim])
//		{
//			CommunicationSumValues(SphericalGravityShellMagneticEnergies[dim], N);
//		}
//	CommunicationSumValues(SphericalGravityShellVolumes, N);

	CommunicationBroadcastValues(SphericalGravityShellMasses, N, ROOT_PROCESSOR);
//Bin Count is only used for debugging.
//	CommunicationSumValues(SphericalGravityShellCellCounts, N);

//	 SphericalGravityMassInterior[0] = SphericalGravityMassShell[0]; //[0] is initialized above.
//	for (int i = 1; i < SphericalGravityActualNumberOfBins; i++)
//	{
////		SphericalGravityMassInterior[i] = SphericalGravityMassInterior[i - 1] + SphericalGravityMassShell[i];
//		SphericalGravityInteriorMasses[i] = SphericalGravityInteriorMasses[i - 1] + SphericalGravityShellMasses[i - 1];
//	}

	arr_delnewset(&SphericalGravityInteriorMasses, N, 0);
	SphericalGravityInteriorMasses[0] = SphericalGravityCentralMass;
	arr_cumsum(SphericalGravityInteriorMasses + 1, SphericalGravityShellMasses, N - 1, SphericalGravityCentralMass);

	switch(SphericalGravityInterpAccelMethod)
	{
	case 1:
		// Calculate the gravity accelerations at the bin edges.
		// SKip the r=0 point.
		for(size_t i = 1; i < N; i++)
			SphericalGravityBinAccels[i] = SphericalGravityConstant * SphericalGravityInteriorMasses[i]
					/ square(SphericalGravityBinLeftEdges[i]);

		break;
	}

//	if(MyProcessorNumber == ROOT_PROCESSOR && MetaData->CycleNumber == 0)
//	{
//		fprintf(stderr, "SphericalGravityActualNumberOfBins = %lld\n", N);
//		arr_printf_pydict("SphericalGravityBinLeftEdges", "%lld:%e", SphericalGravityBinLeftEdges, N);
//		arr_printf_pydict("SphericalGravityShellMasses", "%lld:%e", SphericalGravityShellMasses, N);
//		arr_printf_pydict("SphericalGravityInteriorMasses", "%lld:%e", SphericalGravityInteriorMasses, N);
//		arr_printf_pydict("SphericalGravityBinAccels", "%lld:%e", SphericalGravityBinAccels, N);
//		arr_printf_pydict("SphericalGravityBinAccelSlopes", "%lld:%e", SphericalGravityBinAccelSlopes, N);
//	}

	if(interpolateMissingShells || MetaData->CycleNumber == 0)
	{
		// For shells with no mass (no zones) interpolate the
		// gravity acceleration (g) so it grows smoothly.
		int first_empty, after_first_full;
		FLOAT r1 = 0, r, r2, g1 = 0, g, g2, dgdr;
		for(int j = 0; j < N; j++)
		{
			// Find the next zero shell mass:
			if(SphericalGravityShellMasses[j])
				continue;
			// The first empty shell in a sequence still has
			// the correct g at its left boundary.
			// This is true also for j==0.
			first_empty = j;

			// Find the non-zero shell that follows:
			for(j = first_empty + 1; j < N; j++)
				if(SphericalGravityShellMasses[j])
					break;

			// Make sure we found one and it is not the last one.
			// We don't simply advance j because if the shell to
			// the right hapens to be empty, we are going to miss
			// it on the next loop iteration.
			after_first_full = j + 1;
			if(after_first_full >= N)
				break;

			// Now shell j has g at its left boundary computed
			// with the same enclosed mass as the first_empty
			// shell, which is what we are trying to avoid.
			// Using the g at its right boundary (left of the
			// next shell) to interpolate everything in between.
			g1 = SphericalGravityBinAccels[first_empty];
			r1 = SphericalGravityBinLeftEdges[first_empty];
			g2 = SphericalGravityBinAccels[after_first_full];
			r2 = SphericalGravityBinLeftEdges[after_first_full];
//			TRACEF("  empty=%lld  full=%lld", first_empty, after_first_full);
			dgdr = (g2 - g1) / (r2 - r1);
			for(int k = first_empty + 1; k <= after_first_full; k++)
			{
				r = SphericalGravityBinLeftEdges[k];
				g = g1 + (r - r1) * dgdr;
				SphericalGravityBinAccels[k] = g;
			}
		}
	}

	switch(SphericalGravityInterpAccelMethod)
	{
	case 1:
		// Calculate the slopes of the gravity accelerations between every two bin edges.
		for(size_t i = 0; i < N - 1; i++)
			SphericalGravityBinAccelSlopes[i] = (SphericalGravityBinAccels[i + 1] - SphericalGravityBinAccels[i])
					/ (SphericalGravityBinLeftEdges[i + 1] - SphericalGravityBinLeftEdges[i]);

		SphericalGravityBinAccelSlopes[N - 1] = 0;
		break;
	}

//	if(MyProcessorNumber == ROOT_PROCESSOR && (MetaData->CycleNumber == 0 || debug))
	{
		size_t M = 20;
		TRACEF("SphericalGravity:  GravityConstant = %e    ActualNumberOfBins = %lld    RadiusRange  = %e .. %e",
				SphericalGravityConstant, N, SphericalGravityInnerCutoffRaduis, SphericalGravityOuterCutoffRaduis);
//		arr_printf_pydict("\nSphericalGravityCenter", "%lld:%e", SphericalGravityCenter, 3);
//		arr_printf_pydict("\nSphericalGravityBinLeftEdges", "%lld:%e", SphericalGravityBinLeftEdges, M);
//		arr_printf_pydict("\nSphericalGravityShellMasses", "%lld:%e", SphericalGravityShellMasses, M);
//		arr_printf_pydict("\nSphericalGravityInteriorMasses", "%lld:%e", SphericalGravityInteriorMasses, M);
//		arr_printf_pydict("\nSphericalGravityBinAccels", "%lld:%e", SphericalGravityBinAccels, M);
//		arr_printf_pydict("\nSphericalGravityBinAccelSlopes", "%lld:%e", SphericalGravityBinAccelSlopes, M);
//		fflush(stdout);
	}
//	for (int dim = 0; dim < GridRank; dim++)
//		SphericalGravityCentersOfMass[dim] = arr_sum(SphericalGravityShellCentersOfMass[dim], N);
//
//	SphericalGravityKineticEnergy = 0;
//	for (int dim = 0; dim < GridRank; dim++)
//		SphericalGravityKineticEnergy += SphericalGravityKineticEnergies[dim] = arr_sum(
//				SphericalGravityShellKineticEnergies[dim], N);
//
//	SphericalGravityMagneticEnergy = 0;
//	for (int dim = 0; dim < 3; dim++)
//		if (SphericalGravityShellMagneticEnergies[dim])
//		{
//			SphericalGravityMagneticEnergy +=
//			SphericalGravityMagneticEnergies[dim] = arr_sum(SphericalGravityShellMagneticEnergies[dim], N);
//		}

//	FLOAT burnedR = estimateBurnedRadius();
//	if(burnedR)
//		CommunicationBroadcastValues(&burnedR, 1, MyProcessorNumber);
//	CommunicationR

	return SUCCESS;
}

int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData)
{
	return SphericalGravityComputePotential(LevelArray, MetaData, false);
}

/*
 * Returns the magnitude of the gravity acceleration, g at radius r.
 * Returns -1.0 if error.
 */
float SphericalGravityGetAt(FLOAT r)
{
	if(!UseSphericalGravity || (SphericalGravityInnerCutoffRaduis > 0 && r < SphericalGravityInnerCutoffRaduis)
			|| (SphericalGravityOuterCutoffRaduis >= 0 && r > SphericalGravityOuterCutoffRaduis)
			|| SphericalGravityConstant == 0)
		return 0;

	size_t rbin = SphericalGravityComputeBinIndex(r);
	if(-1 == rbin)
		return 0;

	// SphericalGravityInterpAccelMethod:
	// 0 -- Take the enclosed mass, M_encl, for that bin and return
	//		G * M_encl / r**2
	//      This is always the treatment of the last bin.
	// 1 -- we linearly interpolate g between the edges of the bin.
	//		The outmost bin is considered void of material its enclosed mass
	//		is the entire WD mass. In this case we calculate as in case 0.
	int method = (rbin < SphericalGravityActualNumberOfBins - 1) ? SphericalGravityInterpAccelMethod : 0;
	float g_accel;

	switch(method)
	{
	case 0:
		g_accel = SphericalGravityConstant * SphericalGravityInteriorMasses[rbin] / square(r);
		break;
	case 1:
		// Use linear interpolation between the bin left and right edge.
		// Use the pre-calculated coefficients.
		g_accel = SphericalGravityBinAccels[rbin]
				+ SphericalGravityBinAccelSlopes[rbin] * (r - SphericalGravityBinLeftEdges[rbin]);
		break;
	}

	return g_accel;
}

int SphericalGravityWriteRadialProfile(char* name, TopGridData& MetaData, HierarchyEntry* TopGrid)
{
	HierarchyIterator it = HierarchyIterator(TopGrid);
	for(grid* g = it.firstAtTop(); g; g = it.next())
	{
//		TRACEF(" WriteRadialProfile level %lld", it.currentLevel);
		if(SUCCESS == g->WriteRadialProfile(name, it.currentLevel))
		{
//			return SUCCESS;
		}
	}

	return SUCCESS;
}

int SphericalGravityWritePotential(char * name, TopGridData& MetaData, HierarchyEntry* TopGrid)
{
	SphericalGravityWriteRadialProfile(name, MetaData, TopGrid);

	if(MyProcessorNumber != ROOT_PROCESSOR || UseSphericalGravity == 0 || !SphericalGravityWritePotentialSwitch)
	{
		return SUCCESS;
	}

	if(SphericalGravityInteriorMasses == NULL)
	{
		fprintf(stderr,
				"SphericalGravityWritePotential: No Mass Defined, not writing. This can happen after restart.\n");
		return SUCCESS;
	}

//	typedef short hid_t;
//	hid_t file_id, mass_dset, shell_dset, radius_dset, count_dset, dataspace_id;
//	hid_t ledges_dset, redges_dset, volumes_dset;
//	herr_t status, h5_status, h5_error = -1;
//
//	char filename[100];
//	sprintf(filename, "%s.SphericalGravity.h5", name);
//
//	//file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//	file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
//	if(file_id == -1)
//	{
//		fprintf(stderr, "ERROR IN ERROR: ignore previous warning.  Opening hdf5 file.\n");
//		file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
//	}
//
//	//Create Dataspace
//	hsize_t number[1] = { SphericalGravityActualNumberOfBins };
//	dataspace_id = H5Screate_simple(1, number, NULL);
//
//	//create set
//	//                       above, name,      datatype,  shape of data, Something I dont get
//	mass_dset = H5Dcreate(file_id, "SphericalGravityMassInterior", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	shell_dset = H5Dcreate(file_id, "SphericalGravityMassShell", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	radius_dset = H5Dcreate(file_id, "SphericalGravityRadius", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	ledges_dset = H5Dcreate(file_id, "SphericalGravityLeftEdges", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	redges_dset = H5Dcreate(file_id, "SphericalGravityRightEdges", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	volumes_dset = H5Dcreate(file_id, "SphericalGravityCellVolumes", HDF5_PREC, dataspace_id, H5P_DEFAULT);
//	count_dset = H5Dcreate(file_id, "SphericalGravityBinCount", HDF5_I8, dataspace_id, H5P_DEFAULT);
//
//	//Write the Data Set
//	//                (set, memory type, mem. space, file space, transfer details, actual data)
//	status = H5Dwrite(mass_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityInteriorMasses);
//	status = H5Dwrite(shell_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellMasses);
//	status = H5Dwrite(radius_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinCenters);
//	status = H5Dwrite(ledges_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinLeftEdges);
//	status = H5Dwrite(redges_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinRightEdges);
//	status = H5Dwrite(volumes_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellVolumes);
//	status = H5Dwrite(count_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellCellCounts);
//
//	status = H5Sclose(dataspace_id);
//	status = H5Dclose(mass_dset);
//	status = H5Dclose(shell_dset);
//	status = H5Dclose(radius_dset);
//	status = H5Dclose(ledges_dset);
//	status = H5Dclose(redges_dset);
//	status = H5Dclose(volumes_dset);
//	status = H5Dclose(count_dset);

//	typedef short hid_t;
	hid_t file_id, dataset, dataspace_id;
	hsize_t dims[1];
	herr_t h5_status;

	char filename[FILENAME_MAX];
	sprintf(filename, "%s.sphericalGravity.h5", name);

	if((file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) == -1)
		ENZO_VFAIL("SphericalGravityWritePotential: Error opening hdf5 file `%s'\n", filename)

// A HDF5 writing pattern:
// dataspace_id = H5Screate_simple(rank, dims[], max_dims[])
// dataset = H5Dcreate(output file handle, data set name, data buffer datatype, shape of data handle, H5P_DEFAULT)
// h5_status = H5Dwrite(dset, memory type, buf space, file space, H5P_DEFAULT, data buffer)

// All written arrays have the same dims.
	dims[0] = SphericalGravityActualNumberOfBins;
	dataspace_id = H5Screate_simple(1, dims, NULL);

	dataset = H5Dcreate(file_id, "SphericalGravityLeftEdges", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	h5_status = H5Dwrite(dataset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinLeftEdges);
	h5_status = H5Dclose(dataset);

	dataset = H5Dcreate(file_id, "SphericalGravityMassInterior", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	h5_status = H5Dwrite(dataset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityInteriorMasses);
	h5_status = H5Dclose(dataset);

	if(SphericalGravityInterpAccelMethod == 1)
	{
		dataset = H5Dcreate(file_id, "SphericalGravityBinAccels", HDF5_PREC, dataspace_id, H5P_DEFAULT);
		h5_status = H5Dwrite(dataset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinAccels);
		h5_status = H5Dclose(dataset);

		dataset = H5Dcreate(file_id, "SphericalGravityBinAccelSlopes", HDF5_PREC, dataspace_id, H5P_DEFAULT);
		h5_status = H5Dwrite(dataset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinAccelSlopes);
		h5_status = H5Dclose(dataset);
	}

	h5_status = H5Sclose(dataspace_id);
	h5_status = H5Fclose(file_id);

	return SUCCESS;
}

