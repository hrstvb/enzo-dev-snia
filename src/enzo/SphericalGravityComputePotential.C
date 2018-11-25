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

int CommunicationBroadcastValues(FLOAT *Values, int Number, int BroadcastProcessor);

int SphericalGravityAllocateBins(Eint64** countBins, FLOAT** massBins, FLOAT** centersOfMassBins, FLOAT** kineticEBins,
FLOAT** magneticEBins,
									FLOAT** volumeBins, int domainRank)
{
	size_t &N = SphericalGravityActualNumberOfBins;
	if (countBins)
		arr_delnewset(countBins, N, 0);
	if (massBins)
		for (int dim = 0; dim < domainRank; dim++)
			arr_delnewset(massBins + dim, N, 0);
	if (kineticEBins)
		for (int dim = 0; dim < domainRank; dim++)
			arr_delnewset(kineticEBins + dim, N, 0);
	if (magneticEBins && (UseMHD || UseMHDCT))
		for (int dim = 0; dim < 3; dim++)
			arr_delnewset(magneticEBins + dim, N, 0);
	if (massBins)
		arr_delnewset(massBins, N, 0);
	if (volumeBins)
		arr_delnewset(volumeBins, N, 0);
}

int SphericalGravityDetermineUniformBins()
{
	int nOk = SphericalGravityInnerRadius >= 0;
	nOk += SphericalGravityNumberOfBins > 0;
	nOk += SphericalGravityBinSize > 0;
	nOk += SphericalGravityOuterRadius > 0;
	if (nOk < 3)
	{
		fprintf(stderr, "FAILURE: At least 3 of the following conditions need to be met for uniform bins:"
				"SphericalGravityInnerRadius >= 0, SphericalGravityNumberOfBins > 0, "
				"SphericalGravityBinSize > 0,SphericalGravityOuterRadius > 0,"
				"but the values are %"ISYM", %"FSYM", %"FSYM", %"FSYM".\n",
				SphericalGravityInnerRadius, SphericalGravityNumberOfBins, SphericalGravityBinSize,
				SphericalGravityOuterRadius);
		ENZO_FAIL("SphericalGravityInitializeUniformBins: "
					"Improper parameters inner/outer radius, bin number, bin size.\n");
	}

	if (nOk == 4)
	{
		printf("SphericalGravityInitializeUniformBins: Warning: "
				"Ignored parameter SphericalGravityOuterRadius. "
				"Uniform bis caalculated based on the inner radius, "
				"the bin number and size.\n");
	}

	if (nOk == 4 || SphericalGravityOuterRadius < 0)
	{
		//nothing to do.
	}
	else if (SphericalGravityInnerRadius < 0)
	{
		SphericalGravityInnerRadius = SphericalGravityOuterRadius
				- SphericalGravityNumberOfBins * SphericalGravityBinSize;
		if (SphericalGravityInnerRadius < 0)
			ENZO_FAIL("SphericalGravityInitializeUniformBins: "
						"Improper values fo outer radius, bin number, bin size."
						"The calculated inner radius is negative.\n");
	}
	else
	{
		if (SphericalGravityInnerRadius >= SphericalGravityOuterRadius)
			ENZO_FAIL("SphericalGravityInitializeUniformBins: "
						"Improper parameters inner and outer radii."
						"The inner radius should be < the outer radius.\n");

		if (SphericalGravityBinSize <= 0)
		{
			SphericalGravityBinSize = (SphericalGravityOuterRadius - SphericalGravityInnerRadius)
					/ SphericalGravityNumberOfBins;
		}
		if (SphericalGravityNumberOfBins <= 0)
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

	arr_delnewset(&SphericalGravityBinLeftEdges, N, SphericalGravityBinSize);

	SphericalGravityBinLeftEdges[0] = 0.0;
	if (SphericalGravityHasCentralBin)
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
	if (SphericalGravityUniformBins)
		return SphericalGravityDetermineUniformBins();
//	Initialize non-uniform bins here. This may require communication.
	ENZO_FAIL("Non-uniform spherical gravity bins are not implemented yet.\n");
	return FAIL; //Not implemented
}

size_t SphericalGravityComputeUniformBinIndex(FLOAT r)
{
	if (r < SphericalGravityInnerRadius)
		return 0;
	if (r > SphericalGravityOuterRadius)
		return SphericalGravityActualNumberOfBins - 1;
	size_t i = (r - SphericalGravityInnerRadius) / SphericalGravityBinSize + SphericalGravityHasCentralBin;
	return i;
}

size_t SphericalGravityComputeBinIndex(FLOAT r)
{
	if (SphericalGravityBinSize > 0)
		return SphericalGravityComputeUniformBinIndex(r);
	size_t i = findmaxlte(SphericalGravityBinLeftEdges, SphericalGravityActualNumberOfBins, r);
	return i;
}

int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData)
{
	if (UseSphericalGravity == 0)
		return SUCCESS;

	if (SphericalGravityDetermineBins() == FAIL)
		ENZO_FAIL("Coudn't determine spherical gravity bins.\n");
	SphericalGravityAllocateBins(&SphericalGravityShellCellCounts, SphericalGravityShellCentersOfMass,
									SphericalGravityShellKineticEnergies, SphericalGravityShellMagneticEnergies,
									&SphericalGravityShellMasses, &SphericalGravityShellVolumes, MetaData->TopGridRank);

	LevelHierarchyEntry *Temp;
	Temp = LevelArray[0];
	while (Temp != NULL)
	{
		Temp->GridData->SphericalGravityAddMassToShell();
		Temp = Temp->NextGridThisLevel;
	}

	const size_t &N = SphericalGravityActualNumberOfBins;

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
	arr_cumsum(SphericalGravityShellMasses, SphericalGravityInteriorMasses + 1, N - 1, SphericalGravityCentralMass);

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

	return SUCCESS;
}

int SphericalGravityWritePotential(char * name)
{
	if (MyProcessorNumber != ROOT_PROCESSOR || UseSphericalGravity == 0 || !SphericalGravityWritePotentialSwitch)
	{
		return SUCCESS;
	}

	if (SphericalGravityInteriorMasses == NULL)
	{
		fprintf(stderr, "SphericalGravityWritePotential: No Mass Defined, not writing.\n");
		return SUCCESS;
	}

	typedef short hid_t;
	hid_t file_id, mass_dset, shell_dset, radius_dset, count_dset, dataspace_id;
	hid_t ledges_dset, redges_dset, volumes_dset;
	herr_t status, h5_status, h5_error = -1;

	char filename[100];
	sprintf(filename, "%s.SphericalGravity.h5", name);

	//file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id == -1)
	{
		fprintf(stderr, "ERROR IN ERROR: ignore previous warning.  Opening hdf5 file.\n");
		file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	}

	//Create Dataspace
	hsize_t number[1] = { SphericalGravityActualNumberOfBins };
	dataspace_id = H5Screate_simple(1, number, NULL);

	//create set
	//                       above, name,      datatype,  shape of data, Something I dont get
	mass_dset = H5Dcreate(file_id, "SphericalGravityMassInterior", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	shell_dset = H5Dcreate(file_id, "SphericalGravityMassShell", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	radius_dset = H5Dcreate(file_id, "SphericalGravityRadius", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	ledges_dset = H5Dcreate(file_id, "SphericalGravityLeftEdges", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	redges_dset = H5Dcreate(file_id, "SphericalGravityRightEdges", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	volumes_dset = H5Dcreate(file_id, "SphericalGravityCellVolumes", HDF5_PREC, dataspace_id, H5P_DEFAULT);
	count_dset = H5Dcreate(file_id, "SphericalGravityBinCount", HDF5_I8, dataspace_id, H5P_DEFAULT);

	//Write the Data Set
	//                (set, memory type, mem. space, file space, transfer details, actual data)
	status = H5Dwrite(mass_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityInteriorMasses);
	status = H5Dwrite(shell_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellMasses);
	status = H5Dwrite(radius_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinCenters);
	status = H5Dwrite(ledges_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinLeftEdges);
	status = H5Dwrite(redges_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinRightEdges);
	status = H5Dwrite(volumes_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellVolumes);
	status = H5Dwrite(count_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityShellCellCounts);

	status = H5Sclose(dataspace_id);
	status = H5Dclose(mass_dset);
	status = H5Dclose(shell_dset);
	status = H5Dclose(radius_dset);
	status = H5Dclose(ledges_dset);
	status = H5Dclose(redges_dset);
	status = H5Dclose(volumes_dset);
	status = H5Dclose(count_dset);
	status = H5Fclose(file_id);

	return SUCCESS;
}
