#include "Grid.h"

template<typename T>
bool grid::PointInGridActiveNB(T* point)
{
	for(int dim = 0; dim < GridRank; dim++)
	{
		T p = point[dim];
		if(p <= GridLeftEdge[dim] + CellLeftEdge[dim][NumberOfGhostZones])
			return false;
		if(p >= GridRightEdge[dim] - CellLeftEdge[dim][GridDimension[dim] - NumberOfGhostZones])
			return false;
	}
	return true;
}

bool grid::PointInChildrenActiveNB(FLOAT* point, HierarchyEntry* firstChild)
{
	HierarchyEntry* child = firstChild;
	while(child)
	{
		grid* g = child->GridData;
		if(g && g->PointInGridActiveNB(point))
			return true;
		child = child->NextGridThisLevel;
	}
	return false;
}

bool grid::PointInChildrenActiveNB(FLOAT* point, LevelHierarchyEntry* myLevelHierarchyEntry)
{
	HierarchyEntry* he = (myLevelHierarchyEntry) ? myLevelHierarchyEntry->GridHierarchyEntry : NULL;
	return PointInChildrenActiveNB(point, (he) ? he->NextGridNextLevel : NULL);
}
