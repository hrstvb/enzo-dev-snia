#include "myenzoutils.h"
#include "Grid.h"

bool grid::ijkInGrid(int ijk[])
{
	for(int dim = 0; dim < GridRank; dim++)
	{
		if(ijk[dim] < 0)
			return false;
		if(ijk[dim] >= GridDimension[dim])
			return false;
	}
	return true;
}

bool grid::ijkInGrid(int i, int j, int k)
{
	int ijk[3];// = { i, j, k };
	ijk[0] = i;
	ijk[1] = j;
	ijk[2] = k;
	return ijkInGrid(ijk);
}

template<typename T>
bool grid::PointInGridWithGhost(T* point)
{
	for(int dim = 0; dim < GridRank; dim++)
	{
		T p = point[dim];
		if(p <= CellLeftEdge[dim][0])
			return false;
		if(p >= CellLeftEdge[dim][GridDimension[dim]])
			return false;
	}
	return true;
}

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

void grid::getGhostEdges(FLOAT ledge[], FLOAT redge[])
{
	arr_set(ledge, MAX_DIMENSION, 0);
	arr_set(redge, MAX_DIMENSION, 0);
	for(int dim = 0; dim < GridRank; dim++)
	{
		ledge[dim] = CellLeftEdge[dim][0];
		redge[dim] = CellLeftEdge[dim][GridDimension[dim]];
	}
}

int grid::intersect(FLOAT ledge[], FLOAT redge[])
{
	FLOAT ldomain[MAX_DIMENSION];
	FLOAT rdomain[MAX_DIMENSION];
	getGhostEdges(ldomain, rdomain);
	int retval = intersectRectangles(ledge, redge, ldomain, rdomain, GridRank);
	return retval;
}

int grid::intersectActive(FLOAT ledge[], FLOAT redge[])
{
	int retval = intersectRectangles(ledge, redge, GridLeftEdge, GridRightEdge, GridRank);
	return retval;
}

int grid::intersect(long lijk[], long rijk[])
{
	long long LR[2][MAX_DIMENSION];
	arr_set(LR[0], 2 * MAX_DIMENSION, 0);
	arr_set(LR[1], GridRank, -1);
	arr_xpy(LR[1], GridDimension, GridRank);
	return intersectRectangles(lijk, rijk, LR[0], LR[1], GridRank);
}

int grid::intersectActive(long lijk[], long rijk[])
{
	return intersectRectangles(lijk, rijk, GridStartIndex, GridEndIndex, GridRank);
}
