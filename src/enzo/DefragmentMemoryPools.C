/***********************************************************************
/
/  DEFRAGMENT MEMORYPOOLS
/
/  written by: Tom Abel
/  date:       September 2011
/  modified1:
/
/  PURPOSE:
/   Copy the entire hierarchy from into newly allocated MemoryPools
/   defragmenting the memory used by grids, particles, baryon fields
/   and fluxes. 
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
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
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "MemoryPool.h"


int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);


int DefragmentMemoryPools(LevelHierarchyEntry *LevelArray[])
{

#ifdef MEMORY_POOL
  LCAPERF_START("Defragment Memory Pools"); // 
  if (debug) 
    fprintf(stdout, "Defragmenting memory pools ?.\n");

  // allocate new grid memory pool we will copy to
  int level, i, grid1;
  HierarchyEntry **Grids;
  int NumberOfGrids;
#ifdef HARDPART
  // This will hold the new Level Array
  LevelHierarchyEntry *NewLevelArray[MAX_DEPTH_OF_HIERARCHY];


  int GridObjectSize = sizeof(grid);
  size_t gsize = GridMemoryPool->ReturnUsedMemoryPoolSize();
  size_t ngsize = (ceil(float(gsize)/GridObjectMemorySize)+100);
  Mpool *OldGridPool = NULL;
  OldGridPool = GridObjectMemoryPool;  // copy old one
// create new space and set it to the global pointer 
// this allows the overladed new and delete operators to do the right thing.
  GridObjectMemoryPool = new MPool::MemoryPool(1, GridObjectMemorySize*ngsize, 
					       GridObjectSize,
					       GridObjectMemorySize*ngsize);
    
#endif
#define BCHUNKNUMBER 64
  // do nothing if fragmentation is not bad
  // if (   (BaryonFieldMemoryPool->ReturnFreeMemoryPoolSize())  
  // 	 < (BaryonFieldMemoryPool->ReturnUsedMemoryPoolSize()))
  //   return SUCCESS;

  if (debug) 
    fprintf(stdout, "Defragmenting baryon memory \n");

  size_t bsize = BaryonFieldMemoryPool->ReturnUsedMemoryPoolSize();
  size_t bobjs = sizeof(float);
  size_t nbsize = ceil(float(bsize)*1.5);
  if (debug) 
    fprintf(stdout, "Total: %0.3f Mb, Used: %0.3f Mb, Free: %0.3f Mb = New: %0.3f Mb.\n", 
	    BaryonFieldMemoryPool->ReturnTotalMemoryPoolSize()/1048576.0,
	    BaryonFieldMemoryPool->ReturnUsedMemoryPoolSize()/1048576.0,
	    BaryonFieldMemoryPool->ReturnFreeMemoryPoolSize()/1048576.0,
	    nbsize/1048576.0);

  //  return SUCCESS;

#if 0
  MPool::MemoryPool *OldBaryonPool;
  OldBaryonPool = BaryonFieldMemoryPool;
  BaryonFieldMemoryPool = new MPool::MemoryPool(3, nbsize,
					sizeof(float)*BCHUNKNUMBER,
					nbsize);
  if (debug)
    fprintf(stdout, "allocated new baryon memory pool.\n");
#endif
  
  // Now reallocate all fields. 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->DefragmentBaryonMemoryPool();
      Grids[grid1]->GridData->DefragmentParticleMemoryPool();
    }  // loop for grid1
  }  // loop for level

  //  delete OldBaryonPool; // this will clean up all the old memory pool

#ifdef MPI
// make sure defragmentation is finished for everyone before proceeding
  CommunicationBarrier(); 
#endif

#endif
  return SUCCESS;

}







