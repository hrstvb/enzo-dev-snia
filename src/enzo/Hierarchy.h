/***********************************************************************
/
/  HIERARCHY STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifndef HIERARCHY_ENTRY_DEFINED__
#define HIERARCHY_ENTRY_DEFINED__

struct HierarchyEntry
{
#ifdef MEMORY_POOL
void* operator new(size_t object_size)
{
  return HierarchyEntryMemoryPool->GetMemory(object_size);
}
void operator delete(void* object)
{
  HierarchyEntryMemoryPool->FreeMemory(object);
  return;
}
#endif

  HierarchyEntry *NextGridThisLevel; /* pointer to the next grid on level */
  HierarchyEntry *NextGridNextLevel; /* pointer to first child of this grid */
  HierarchyEntry *ParentGrid;        /* pointer to this grid's parent */
  grid           *GridData;          /* pointer to this grid's data */
};



#endif
