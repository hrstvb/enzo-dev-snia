#include "myenzoutils.h"
#include "DebugMacros.h"
#include "ErrorExceptions.h"
#include "LevelArrayIterator.h"
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level);

LevelArrayIterator::LevelArrayIterator(LevelHierarchyEntry** levelArray)
{
	this->levelArray = levelArray;
}

int LevelArrayIterator::countLevels()
{
	for(int n = 0; n < MAX_DEPTH_OF_HIERARCHY; n++)
		if(levelArray[n] == NULL)
			return n;
	return MAX_DEPTH_OF_HIERARCHY;
}

grid* LevelArrayIterator::getCurrentGrid(grid** parent)
{
	if(currentEntry)
	{
		if(parent)
			*parent = getCurrentParent();
		return currentEntry->GridData;
	}

	if(parent)
		*parent = NULL;
	return NULL;;
}

grid* LevelArrayIterator::getCurrentParent()
{
	if(currentEntry)
	{
		HierarchyEntry* he = currentEntry->GridHierarchyEntry;
		if(he)
			if(he = he->ParentGrid)
				return he->GridData;
	}
	return NULL;
}

grid* LevelArrayIterator::firstFromTop()
{
	return firstFromTop(NULL);
}

grid* LevelArrayIterator::firstFromTop(grid** parent)
{
	currentLevel = 0;
	currentEntry = levelArray[0];
	forward = true;
	endOfLevel = false;
	return getCurrentGrid(parent);
}

grid* LevelArrayIterator::firstFromFinest()
{
	return firstFromFinest(NULL);
}

grid* LevelArrayIterator::firstFromFinest(grid** parent)
{
	forward = false;
	endOfLevel = false;
	currentLevel = countLevels() - 1;
	if(currentLevel >= 0)
		currentEntry = levelArray[currentLevel];
	else
		currentEntry = NULL;
	return getCurrentGrid(parent);
}

grid* LevelArrayIterator::firstFromLevel(int level)
{
	return firstFromLevel(level, NULL);
}

grid* LevelArrayIterator::firstFromLevel(int level, grid** parent)
{
	currentLevel = level;
	currentEntry = levelArray[level];
	forward = true;
	endOfLevel = currentEntry == NULL;
	return getCurrentGrid(parent);
}

grid* LevelArrayIterator::next()
{
	return next(NULL);
}

grid* LevelArrayIterator::next(grid** parent)
{
	if(currentEntry)
	{
		currentEntry = currentEntry->NextGridThisLevel;
		if(currentEntry == NULL)
		{
			currentLevel += (forward) ? 1 : -1;
			if(0 <= currentLevel && currentLevel < MAX_DEPTH_OF_HIERARCHY)
				currentEntry = levelArray[currentLevel];
		}
	}

	return getCurrentGrid(parent);
}

grid* LevelArrayIterator::nextThisLevel()
{
	return nextThisLevel(NULL);
}

grid* LevelArrayIterator::nextThisLevel(grid** parent)
{
	LevelHierarchyEntry* lhe = NULL;
	if(!endOfLevel && currentEntry != NULL)
		lhe = currentEntry->NextGridThisLevel;

	if(endOfLevel = (lhe == NULL))
	{
		if(parent)
			*parent == NULL;
		return NULL;
	}

	currentEntry = lhe;
	return getCurrentGrid(parent);
}

grid* LevelArrayIterator::prev()
{
	return prev(NULL);
}

grid* LevelArrayIterator::prev(grid** parent)
{
	if(currentEntry == NULL || currentEntry == levelArray[currentLevel])
	{
		currentLevel += (forward) ? -1 : 1;
		currentEntry = NULL;
	}

	if(currentLevel < 0 || currentLevel >= MAX_DEPTH_OF_HIERARCHY)
	{
		if(parent)
			parent == NULL;
		return NULL;
	}

	LevelHierarchyEntry* lhe = levelArray[currentLevel];
	while(lhe)
	{
		if(lhe->NextGridThisLevel == currentEntry)
			break;
		lhe = lhe->NextGridThisLevel;
	}
	currentEntry = lhe;

	return getCurrentGrid(parent);
}

HierarchyIterator::HierarchyIterator(const HierarchyEntry* const topGrid) :
		topGrid(topGrid)
{
	setCurrent(NULL, -1, NULL);
}

//grid* HierarchyIterator::firstAtTop()
//{
//	if(topGrid == 0)
//		return NULL;
//
//	currentHEntry = topGrid;
//	currentLevel = 0;
//	return topGrid->GridData;
//}
//
//grid* HierarchyIterator::next()
//{
//	if(currentHEntry == NULL)
//		return NULL;
//
//	HierarchyEntry *he2;
//
//	// Try going down
//	he2 = currentHEntry->NextGridNextLevel;
//	if(he2)
//	{
//		currentLevel++;
//		currentHEntry = he2;
//		return he2->GridData;
//	}
//
//	// Try going right, if can't -- go up and try right again.
//	while(currentHEntry)
//	{
//		he2 = currentHEntry->NextGridThisLevel;
//		if(he2)
//		{
//			currentHEntry = he2;
//			return he2->GridData;
//		}
//
//		currentHEntry = currentHEntry->ParentGrid;
//		currentLevel--;
//	}
//
//	return NULL;
//}

grid* HierarchyIterator::getCurrentParent()
{
	if(currentHEntry)
	{
		const HierarchyEntry* const he = currentHEntry->ParentGrid;
		if(he)
			return he->GridData;
	}
	return NULL;
}

grid* HierarchyIterator::setCurrent(const HierarchyEntry* const newEntry, const int newLevel, grid** const parent)
{
	currentHEntry = newEntry;
	currentLevel = newLevel;
	if(parent)
		*parent = getCurrentParent();
	return (currentHEntry) ? currentHEntry->GridData : NULL;
}

grid* HierarchyIterator::first(bool parentsAfterChildren, bool parentsBeforeChildren)
{
	this->parentsAfterChildren = parentsAfterChildren;
	this->parentsBeforeChildren = parentsBeforeChildren;
	if(topGrid == NULL)
		return setCurrent(NULL, -1, NULL);

	currentHEntry = topGrid;
	currentLevel = 0;
	backing = isLeaf = topGrid->NextGridNextLevel == NULL;
	if(parentsBeforeChildren || isLeaf)
		return setCurrent(topGrid, 0, NULL);

	return next();
}

grid* HierarchyIterator::firstAtTop()
{
	return first(false, true);
}

grid* HierarchyIterator::firstAtFinest()
{
	return first(true, false);
}

grid* HierarchyIterator::next()
{
	if(currentHEntry == NULL)
		return setCurrent(NULL, -1, NULL);

	const HierarchyEntry* he = currentHEntry;
	const HierarchyEntry* he2;
	int l2 = currentLevel;

	if(backing)
	{
		while(he)
		{
			he2 = he->NextGridThisLevel;
			if(he2)
			{
				backing = false;
				he = he2;
				break;
			}
			he = he->ParentGrid;
			l2--;
			if(parentsAfterChildren && he)
				return setCurrent(he, l2, NULL);
		}
	}
	else
	{
		he = currentHEntry->NextGridNextLevel;
		l2++;
	}

	if(he == NULL)
		return setCurrent(NULL, -1, NULL);

	while(he)
	{
		he2 = he->NextGridNextLevel;
		isLeaf = he2 == NULL;
		if(isLeaf)
		{
			backing = true;
			return setCurrent(he, l2, NULL);
		}
		else if(parentsBeforeChildren)
			return setCurrent(he, l2, NULL);
		he = he2;
		l2++;
	}

	return setCurrent(NULL, -1, NULL);
}

RebuildHierarchyIterator::RebuildHierarchyIterator(int maxRefineLevel, HierarchyEntry* topGrid, TopGridData* metaData)
{
	this->maxRefLevel = maxRefineLevel;
	maxRefLevel = (maxRefLevel < 0) ? 0 : maxRefLevel;
	maxRefLevel = (maxRefLevel < MAX_DEPTH_OF_HIERARCHY) ? maxRefLevel : MAX_DEPTH_OF_HIERARCHY - 1;

	this->topGrid = topGrid;
	this->metaData = metaData;
	currentEntry = NULL;
	currentLevel = -1;
	startingNewLevel = false;
	arr_set(levelArray, MAX_DEPTH_OF_HIERARCHY, 0);
}

grid* RebuildHierarchyIterator::first()
{
	AddLevel(levelArray, topGrid, 0);
	currentEntry = levelArray[currentLevel = 0];
	startingNewLevel = true;
	return currentEntry->GridData;
}

grid* RebuildHierarchyIterator::next()
{
	startingNewLevel = false;

	if(currentEntry == NULL || currentLevel < 0 || currentLevel > maxRefLevel)
		return NULL;

	LevelHierarchyEntry* lhe = currentEntry->NextGridThisLevel;
	if(lhe)
	{
		currentEntry = lhe;
		return lhe->GridData;
	}

	if(currentLevel == maxRefLevel)
		return NULL;

	int MaximumRefinementLevel_original = MaximumRefinementLevel;
	MaximumRefinementLevel = currentLevel + 1;
	TRACEF("Refining hierarchy level %lld ...", currentLevel);
	if(RebuildHierarchy(metaData, levelArray, currentLevel) == FAIL)
		throw EnzoFatalException("Error in RebuildHierarchy.", __FILE__, __LINE__);
	TRACEF("Refining hierarchy level %lld DONE.", currentLevel);
	MaximumRefinementLevel = MaximumRefinementLevel_original;
	currentEntry = levelArray[++currentLevel];
	startingNewLevel = true;
	return (currentEntry) ? currentEntry->GridData : NULL;
}

SiblingIterator::SiblingIterator(const HierarchyEntry* const firstSibling) :
		firstSibling { firstSibling }, currentSibling { NULL }
{
}

SiblingIterator::SiblingIterator(const LevelHierarchyEntry* const firstSibling) :
		firstSibling { (firstSibling) ? firstSibling->GridHierarchyEntry : NULL }, currentSibling { NULL }
{
}

SiblingIterator SiblingIterator::NewFromParent(const HierarchyEntry* const parent)
{
	return SiblingIterator((parent) ? parent->NextGridNextLevel : NULL);
}

SiblingIterator SiblingIterator::NewFromParent(const LevelHierarchyEntry* const parent)
{
	return NewFromParent((parent) ? parent->GridHierarchyEntry : NULL);
}

grid* SiblingIterator::first()
{
	if((currentSibling = firstSibling) == NULL)
		return NULL;
	grid* g = currentSibling->GridData;
	return (g) ? g : next();
}

grid* SiblingIterator::next()
{
	while(currentSibling)
	{
		currentSibling = currentSibling->NextGridThisLevel;
		grid* g = currentSibling->GridData;
		if(g)
			return g;
	}
	return NULL;
}
