/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Tom Abel
/  date:       October 2011
/  modified:   
/  date:       
/              
/
/  PURPOSE: Allocate the global memory pools.
/     
/
************************************************************************/
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void AllocateMemoryPools() 
{
  /* The initial size of the memory pool in units of photon packages.
     Increase the memory pool by 1/4th of the initial size as more
     memory is needed. */
  
  
#ifdef MEMORY_POOL
#ifdef GRID_MEMORY_POOL
  const int GridObjectMemorySize = MEMORY_POOL_SIZE;
  int GridObjectSize = sizeof(grid);
  GridObjectMemoryPool = new MPool::MemoryPool(1, GridObjectMemorySize*GridObjectSize,
					       GridObjectSize,
					       GridObjectMemorySize*GridObjectSize/4, true);
#endif
#ifdef PROTOSUBGRID_MEMORY_POOL
  const int ProtoSubgridMemorySize = 1024*128;
  int ProtoSubgridSize = sizeof(ProtoSubgrid);
  ProtoSubgridMemoryPool = new MPool::MemoryPool(2, ProtoSubgridMemorySize*ProtoSubgridSize,
						 ProtoSubgridSize,
						 ProtoSubgridMemorySize*ProtoSubgridSize, true);
#endif
#ifdef HIERARCHY_MEMORY_POOL
  const int HierarchyEntryMemorySize =1024*128;
  int HierarchyEntrySize = sizeof(HierarchyEntry);
  HierarchyEntryMemoryPool = new MPool::MemoryPool(3, HierarchyEntryMemorySize*HierarchyEntrySize,
						   HierarchyEntrySize,
						   HierarchyEntryMemorySize*HierarchyEntrySize, true);
#endif
  FlaggingFieldMemoryPool = new MPool::MemoryPool(4, sizeof(int)*1024*1024/NumberOfProcessors,
						  sizeof(int)*32,
						  sizeof(int)*200000, true);
  
  ParticleMemoryPool = new MPool::MemoryPool(5, sizeof(FLOAT)*4000000/NumberOfProcessors,
					     sizeof(FLOAT)*64,
					     sizeof(FLOAT)*2000000, true);
  
  BaryonFieldMemoryPool = new MPool::MemoryPool(6, sizeof(float)*24000000/NumberOfProcessors,
						sizeof(float)*64,
						sizeof(float)*6000000, true);

#ifdef TRANSFER
  const int PhotonMemorySize = MEMORY_POOL_SIZE;
  int PhotonSize = sizeof(PhotonPackageEntry);
  PhotonMemoryPool = new MPool::MemoryPool(4, PhotonMemorySize*PhotonSize,
					   PhotonSize,
					   PhotonMemorySize*PhotonSize/4);
#endif
#endif

  return;
}
