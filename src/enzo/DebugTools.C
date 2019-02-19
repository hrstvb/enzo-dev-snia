#include <string.h>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "hdf5.h"
#include <sys/time.h>
#include <stdlib.h>

#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "CommunicationUtilities.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#include "DebugTools.h"
using namespace std;

void RecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count);

//size_t sprintHierarchy(char* s, LevelHierarchyEntry** levelArray)
//{
//	char* t = s;
//	t += sprintf(t, "BEGIN HIERARCHY (on #%lld) ----------------\n", MyProcessorNumber);
//	for(int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
//	{
//		LevelHierarchyEntry* lhe = levelArray[level];
//		while(lhe)
//		{
//			grid* g = lhe->GridData;
//			int gProc = g->ReturnProcessorNumber();
//			grid* pg = NULL;
//			HierarchyEntry* he = lhe->GridHierarchyEntry->ParentGrid;
//			if(he)
//				pg = he->GridData;
//
//			t += sprintf(t, "HIERARCHY (on #%lld) level %" ISYM "    grid %" ISYM "(%p)", MyProcessorNumber, level,
//							g->GetGridID(), g);
//			if(gProc == MyProcessorNumber)
//			{
//				t += sprintf(t, "(local)");
//			}
//			else
//			{
//				t += sprintf(t, "(on #%lld)", gProc);
//			}
//
//			if(pg)
//			{
//				t += sprintf(t, "    parent %" ISYM "%p", pg->GetGridID(), pg);
//			}
//
//			t += sprintf(t, "\n");
//
//			lhe = lhe->NextGridThisLevel;
//		}
//	}
//
//	t += sprintf(t, "END HIERARCHY (on #%lld) ----------------\n", MyProcessorNumber);
//	return t - s;
//}

int sprintHierarchy(char* s, size_t size, LevelHierarchyEntry** levelArray)
{
	size_t total_length = 0;
	int n;

	snlprintf(s, size, &total_length, "BEGIN HIERARCHY (on #%lld) ----------------\n", MyProcessorNumber);
	for(int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	{
		LevelHierarchyEntry* lhe = levelArray[level];
		while(lhe)
		{
			grid* g = lhe->GridData;
			int gProc = g->ReturnProcessorNumber();
			grid* pg = NULL;
			HierarchyEntry* he = lhe->GridHierarchyEntry->ParentGrid;
			if(he)
				pg = he->GridData;

			n = snlprintf(s, size, &total_length, "HIERARCHY (on #%lld) level %" ISYM "    grid %" ISYM "(%p)",
							MyProcessorNumber, level, g->GetGridID(), g);
			if(n < 0) return n;

			if(gProc == MyProcessorNumber)
			{
				n = snlprintf(s, size, &total_length, "(local)");
				if(n < 0) return n;
			}
			else
			{
				n = snlprintf(s, size, &total_length, "(on #%lld)", gProc);
				if(n < 0) return n;
			}

			if(pg)
			{
				n = snlprintf(s, size, &total_length, "    parent %" ISYM "(%p)", pg->GetGridID(), pg);
				if(n < 0) return n;
			}

			n = snlprintf(s, size, &total_length, "\n");
			if(n < 0) return n;

			lhe = lhe->NextGridThisLevel;
		}
	}

	n = snlprintf(s, size, &total_length, "END HIERARCHY (on #%lld) ----------------\n", MyProcessorNumber);
	if(n < 0) return n;

	return total_length;
}

void printHierarchy(LevelHierarchyEntry** levelArray)
{
	const size_t SIZE = 1024;
	char buf[SIZE];
	char* s = buf;
	int n = sprintHierarchy(s, SIZE, levelArray);
	if(n >= SIZE)
	{
		s = new char[n];
		sprintHierarchy(s, n, levelArray);
	}
	fprintf(stderr, s);
}

void printHierarchy0(LevelHierarchyEntry** levelArray)
{
	if(MyProcessorNumber == ROOT_PROCESSOR)
		printHierarchy(levelArray);
}

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid, TopGridData *MetaData)
{

	if(TracerParticlesAddToRestart == FALSE)
		return SUCCESS;

	//Adds evenly spaced tracer particles to an existing run (most likely
	//to a restart from a previous enzo version or otherwise generated dataset.)
	//Requires the following flags to be set:
	//TracerParticlesAddToRestart = 1 (this is not written out, only a trigger)
	//TracerParticleOn = 1
	//TracerParticleCreationSpacing = 0.125 (or whatever)
	//TracerParticleCreationLeftEdge = 0 0 0
	//TracerParticleCreationRightEdge= 1 1 1

	char line[MAX_LINE_LENGTH];
	int dim;

	// Set default values for parameters

	//  Declared in global_data.h (RH)
	//  FLOAT TracerParticleCreationLeftEdge[MAX_DIMENSION];
	//  FLOAT TracerParticleCreationRightEdge[MAX_DIMENSION];
	//  FLOAT TracerParticleCreationSpacing = FLOAT_UNDEFINED;

	TracerParticleCreationSpacing = FLOAT_UNDEFINED;

	for(dim = 0; dim < MAX_DIMENSION; dim++)
	{
		TracerParticleCreationLeftEdge[dim] = FLOAT_UNDEFINED;
		TracerParticleCreationRightEdge[dim] = FLOAT_UNDEFINED;
	}

	// Read until out of lines

	FILE * fptr = fopen(filename, "r");

	while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
	{

		// Read tracer particle parameters

		sscanf(line, "TracerParticleCreationSpacing = %"PSYM, &TracerParticleCreationSpacing);
		sscanf(line, "TracerParticleCreationLeftEdge = %"PSYM" %"PSYM" %"PSYM, TracerParticleCreationLeftEdge,
				TracerParticleCreationLeftEdge + 1, TracerParticleCreationLeftEdge + 2);
		sscanf(line, "TracerParticleCreationRightEdge = %"PSYM" %"PSYM" %"PSYM, TracerParticleCreationRightEdge,
				TracerParticleCreationRightEdge + 1, TracerParticleCreationRightEdge + 2);

	}
	fclose(fptr);

	if(debug)
	{
		fprintf(stderr, "TracerParticleCreation = %"ISYM"\n", MetaData->CycleNumber);
		fprintf(stderr, "TracerParticleCreationSpacing = %"PSYM"\n", TracerParticleCreationSpacing);
		fprintf(stderr, "TracerParticleCreationLeftEdge = %"PSYM" %"PSYM" %"PSYM"\n", TracerParticleCreationLeftEdge[0],
				TracerParticleCreationLeftEdge[1], TracerParticleCreationLeftEdge[2]);
		fprintf(stderr, "TracerParticleCreationRightEdge = %"PSYM" %"PSYM" %"PSYM"\n",
				TracerParticleCreationRightEdge[0], TracerParticleCreationRightEdge[1],
				TracerParticleCreationRightEdge[2]);
	}

	HierarchyEntry * Temp = TopGrid;
	if(TracerParticleCreationSpacing > 0)
	{
		while(Temp != NULL)
		{
			fprintf(stderr, "adding tracer particles\n");
			if(Temp->GridData->TracerParticleCreateParticles(TracerParticleCreationLeftEdge,
																TracerParticleCreationRightEdge,
																TracerParticleCreationSpacing,
																MetaData->NumberOfParticles) == FAIL)
			{
				fprintf(stderr, "Error in grid->TracerParticleCreateParticles.\n");
				ENZO_FAIL("");
			}
			Temp = Temp->NextGridThisLevel;
		}
	}
	// Get the global particle count

	int LocalNumberOfParticles;

	Temp = TopGrid;
	while(Temp != NULL)
	{

		LocalNumberOfParticles = 0;
		LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
		if(debug)
			printf("OldLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles);

#ifdef USE_MPI
		CommunicationAllReduceValues(&LocalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */
		Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);

		LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
		if(debug)
			printf("NewLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles);

		Temp = Temp->NextGridThisLevel;
	}

	Temp = TopGrid;
	PINT ParticleCount = 0;
	RecursivelySetParticleCount(Temp, &ParticleCount);
	MetaData->NumberOfParticles = ParticleCount;

	return SUCCESS;
}

//Writes an HDF5 cube.  Quick and dirty.
// WriteCube(array, [nx,ny,nz], "ID string", dNum, gNum)
// prints out ID string to std err, along with file name.
// Filename is data111(dNum).grid(gNum)
// datasets are all called "monkey"
// Flow:
// 1.) create filename, print message
// 2.) define size of float
// 3.) create file
// 3.5) invert array dimensions
// 4.) create dataspace, set
// 5.) close file,space,set.
void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label)
{

	hid_t file_id, dataset_id, dataspace_id, float_type_id;
	herr_t status, h5_status, h5_error = -1;
	int FieldRankOut = 3;
	hsize_t DimsInv[FieldRankOut];

	char filename[20];

	sprintf(filename, "data111%4.4d.grid%4.4d", dNum, gNum);
	fprintf(stderr, "GPFS WriteCube: %s %s [%"ISYM",%"ISYM",%"ISYM"]\n", string, filename, Dims[0], Dims[1], Dims[2]);

#define floatdcc double
	int jj = sizeof(floatdcc);
	switch(jj)
	{
	case 0:
		float_type_id = H5T_NATIVE_INT;
		break;
	case 4:
		float_type_id = HDF5_R4;
		break;
	case 8:
		float_type_id = HDF5_R8;
		break;
	case 16:
		float_type_id = HDF5_R16;
		break;
	default:
		printf("INCORRECT FLOAT DEFINITION\n");
	}

	//file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if(file_id == -1)
	{
		fprintf(stderr, "ERROR IN ERROR: ignore previous warning.  Opening hdf5 file.\n");
		file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	}
	for(int moo = 0; moo < FieldRankOut; moo++)
		DimsInv[moo] = Dims[FieldRankOut - moo - 1];

	//Create Dataspace
	dataspace_id = H5Screate_simple(FieldRankOut, DimsInv, NULL);

	//create set
	//                       duh, name,      datatype,  shape of data, Something I dont get
	dataset_id = H5Dcreate(file_id, label, float_type_id, dataspace_id, H5P_DEFAULT);

	//Write the Data Set
	//                (set, memory type, mem. space, file space, transfer details, actual data)
	fprintf(stderr, "Writing set\n");
	status = H5Dwrite(dataset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

	status = H5Sclose(dataspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);

}
