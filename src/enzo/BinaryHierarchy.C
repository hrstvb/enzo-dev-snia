/***********************************************************************
/
/  BINARY HIERARCHY STRUCTURE HELPERS
/
/  written by: Matthew Turk
/  date:       October, 2009
/  modified1:
/
/  PURPOSE: Initialize storage inside a hierarchy array container
/
/  REQUIRES: macros_and_parameters.h
/
************************************************************************/

#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "h5utilities.h"
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "BinaryHierarchy.h"

hierarchy_arrays HierarchyArrays;
FLOAT *AllocateNewParticleMemory(int size);
int *AllocateNewFlaggingField(int size);
void FreeParticleMemory(void *PF);
void FreeFlaggingFieldMemory(int *FF);


void InitializeHierarchyArrayStorage(int grid_count)
{
    HierarchyArrays.current_parent = -1;
    HierarchyArrays.grid_count = grid_count;
    HierarchyArrays.ActiveDimensions = AllocateNewFlaggingField(MAX_DIMENSION * grid_count);
    HierarchyArrays.LeftEdges = AllocateNewParticleMemory(MAX_DIMENSION * grid_count);
    HierarchyArrays.RightEdges = AllocateNewParticleMemory(MAX_DIMENSION * grid_count);
    HierarchyArrays.Level = AllocateNewFlaggingField(grid_count);
    HierarchyArrays.ParentIDs = AllocateNewFlaggingField(grid_count);
    HierarchyArrays.Processor = AllocateNewFlaggingField(grid_count);
    HierarchyArrays.NumberOfParticles = AllocateNewFlaggingField(grid_count);
}

void WriteHierarchyArrayStorage(const char* name)
{
    hid_t file_id, group_id, dset_id, file_dsp_id, mem_dsp_id;
    hsize_t dims[2]; /* At most two */
    herr_t h5_error = -1;
    Eint32 ndims;

    if(MyProcessorNumber != ROOT_PROCESSOR) return;
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id == h5_error) ENZO_VFAIL("Failure creating %s", name)

    /* We have the file, now we manually write out each data set */

    /* ActiveDimensions */
    ndims = 2; dims[0] = HierarchyArrays.grid_count; dims[1] = MAX_DIMENSION;
    writeArrayDataset(file_id, HDF5_INT, ndims, dims, "ActiveDimensions", HierarchyArrays.ActiveDimensions);
    writeArrayDataset(file_id, HDF5_PREC, ndims, dims, "LeftEdges", HierarchyArrays.LeftEdges);
    writeArrayDataset(file_id, HDF5_PREC, ndims, dims, "RightEdges", HierarchyArrays.RightEdges);
    
    ndims = 1;
    writeArrayDataset(file_id, HDF5_INT, ndims, dims, "Level", HierarchyArrays.Level);
    writeArrayDataset(file_id, HDF5_INT, ndims, dims, "ParentIDs", HierarchyArrays.ParentIDs);
    writeArrayDataset(file_id, HDF5_INT, ndims, dims, "Processor", HierarchyArrays.Processor);
    writeArrayDataset(file_id, HDF5_INT, ndims, dims, "NumberOfParticles", HierarchyArrays.NumberOfParticles);

    H5Fclose(file_id);
}

void FinalizeHierarchyArrayStorage()
{
  FreeFlaggingFieldMemory(HierarchyArrays.ActiveDimensions);
  FreeParticleMemory((void*) HierarchyArrays.LeftEdges);
  FreeParticleMemory((void*) HierarchyArrays.RightEdges);
  FreeFlaggingFieldMemory(HierarchyArrays.Level);
  FreeFlaggingFieldMemory(HierarchyArrays.ParentIDs);
  FreeFlaggingFieldMemory(HierarchyArrays.Processor);
  FreeFlaggingFieldMemory(HierarchyArrays.NumberOfParticles);
}
