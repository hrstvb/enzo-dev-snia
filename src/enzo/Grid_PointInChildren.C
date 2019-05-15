#include "myenzoutils.h"
#include "Grid.h"

template<typename T>
bool grid::PointInGridActiveNB(T* point)
{
	for(int dim = 0; dim < GridRank; dim++)
	{
		T p = point[dim];
		if(p <= GridLeftEdge[dim])
			return false;
		if(p >= GridRightEdge[dim])
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

void grid::getDomainEdges(FLOAT ledge[], FLOAT redge[])
{
	arr_set(ledge, MAX_DIMENSION, 0);
	arr_set(redge, MAX_DIMENSION, 0);
	for(int dim = 0; dim < GridRank; dim++)
	{
		ledge[dim] = CellLeftEdge[dim][0];
		redge[dim] = CellLeftEdge[dim][GridDimension[dim]];
	}
}

bool grid::intersectDomain(FLOAT ledge[], FLOAT redge[])
{
	FLOAT ldomain[MAX_DIMENSION];
	FLOAT rdomain[MAX_DIMENSION];
	getDomainEdges(ldomain, rdomain);
	return intersectRectangles(ledge, redge, ldomain, rdomain, GridRank);
}
