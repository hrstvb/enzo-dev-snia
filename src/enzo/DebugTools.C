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
#include "LevelArrayIterator.h"
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

int snlprintHierarchy(char* s, size_t size, size_t* length, LevelHierarchyEntry** levelArray,
	const char* const filename, const int linenum)
{
	size_t l0 = *length;
	int n;

	if(filename)
		n = snlprintf(s, size, length, "/- BEGIN HIERARCHY (%s:%d#%lld) ----------------\n", filename, linenum,
						MyProcessorNumber);
	else
		n = snlprintf(s, size, length, "/- BEGIN HIERARCHY (#%lld) ----------------\n", MyProcessorNumber);
	if(n < 0)
		return n;

	LevelArrayIterator it = LevelArrayIterator(levelArray);
//	for(int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
//	{
//		LevelHierarchyEntry* lhe = levelArray[level];
//		while(lhe)
	grid* pg = NULL;
	for(grid* g = it.firstFromTop(&pg); g; g = it.next(&pg))
	{
		int gProc = g->ReturnProcessorNumber();
		int dx2pdx[MAX_DIMENSION];
		if(pg)
			pg->ComputeRefinementFactors(g, dx2pdx);
		else
			arr_set(dx2pdx, MAX_DIMENSION, 1);

		n = snlprintf(s, size, length, "|  HIERARCHY (on #%lld) level %" ISYM
		"    grid %" ISYM "(%p)",
						MyProcessorNumber, it.currentLevel, g->GetGridID(), g);
		if(n < 0)
			return n;

		if(gProc == MyProcessorNumber)
		{
			n = snlprintf(s, size, length, "(local)");
			if(n < 0)
				return n;
		}
		else
		{
			n = snlprintf(s, size, length, "(on #%lld)", gProc);
			if(n < 0)
				return n;
		}

		if(pg)
		{
			n = snlprintf(s, size, length, "   x%lld    parent %" ISYM "(%p)", dx2pdx[0], pg->GetGridID(), pg);
			if(n < 0)
				return n;
		}

		n = snlprintf(s, size, length, "\n");
		if(n < 0)
			return n;
	}
//	}

	if(filename)
		n = snlprintf(s, size, length, "\\- END HIERARCHY (%s:%d#%lld) ----------------\n", filename, linenum,
						MyProcessorNumber);
	else
		n = snlprintf(s, size, length, "\\- END HIERARCHY (#%lld) ----------------\n", MyProcessorNumber);
	if(n < 0)
		return n;

	return *length - l0;
}

int printHierarchy(LevelHierarchyEntry** levelArray, const char* const filename, const int linenum)
{
	const size_t SIZE = 1024;
	char S[SIZE];
	size_t len = 0;
	int n = snlprintHierarchy(S, SIZE, &len, levelArray, filename, linenum);
	if(n < 0)
		return n;

	size_t size = len + 1;
	if(size <= SIZE)
		return fprintf(stderr, S);

	char* s = new char[size];
	len = 0;
	n = snlprintHierarchy(s, size, &len, levelArray, filename, linenum);
	if(n >= 0)
		n = fprintf(stderr, s);
	delete s;
	return n;
}

int printHierarchy0(LevelHierarchyEntry** levelArray, const char* const filename, const int linenum)
{
	if(MyProcessorNumber != ROOT_PROCESSOR)
		return 0;
	return printHierarchy(levelArray, filename, linenum);
}

int printHierarchy(LevelHierarchyEntry** levelArray)
{
	return printHierarchy(levelArray, NULL, -1);
}

int printHierarchy0(LevelHierarchyEntry** levelArray)
{
	return printHierarchy0(levelArray, NULL, -1);
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
