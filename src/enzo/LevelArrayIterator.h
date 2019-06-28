/*
 * LevelArrayIterator.h
 *
 *  Created on: Jan 23, 2019
 *      Author: cdev
 */

#ifndef SRC_ENZO_LEVELARRAYITERATOR_H_
#define SRC_ENZO_LEVELARRAYITERATOR_H_

class grid;
struct TopGridData;
#include "global_data.h"
#include "LevelHierarchy.h"

#ifndef MAX_DEPTH_OF_HIERARCHY
#define MAX_DEPTH_OF_HIERARCHY (50)
#endif

struct LevelArrayIterator
{
	LevelHierarchyEntry** levelArray; //
	int nLevels = -1; ///
	bool forward = true; //
	int currentLevel = 0; //
	LevelHierarchyEntry* currentEntry = NULL;
	bool endOfLevel = false;

	LevelArrayIterator(LevelHierarchyEntry** levelArray); //
	int countLevels();
	grid* getCurrentGrid(grid** parent);
	grid* getCurrentParent();
	grid* firstFromTop();
	grid* firstFromTop(grid** parent);
	grid* firstFromLevel(int level);
	grid* firstFromLevel(int level, grid** parent);
	grid* firstFromFinest();
	grid* firstFromFinest(grid** parent);
	grid* next();
	grid* next(grid** parent);
	grid* nextThisLevel();
	grid* nextThisLevel(grid** parent);
	grid* prev();
	grid* prev(grid** parent);
	int projectChildrenToParents(bool projectB, bool projectE);
};

struct HierarchyIterator
{
	const HierarchyEntry* const topGrid;
	const HierarchyEntry* currentHEntry;int currentLevel;
	bool backing = false;
	bool parentsAfterChildren = true;
	bool parentsBeforeChildren = true;
	bool isLeaf = false;

	HierarchyIterator(const HierarchyEntry* const topGrid);
	grid* first(bool parentsAfterChildren, bool parentsBeforeChildren);
	grid* firstAtTop();
	grid* firstAtFinest();
	grid* next();
	grid* getCurrentParent();
	grid* setCurrent(const HierarchyEntry* const newCurrentEntry, const int newCurrentLevel, grid** const parent);
};

struct RebuildHierarchyIterator
{
	int maxRefLevel;
	HierarchyEntry* topGrid;
	TopGridData* metaData; //
	int currentLevel;
	LevelHierarchyEntry* currentEntry; //
	LevelHierarchyEntry* levelArray[MAX_DEPTH_OF_HIERARCHY];
	bool startingNewLevel;

	RebuildHierarchyIterator(int maxRefineLevel, HierarchyEntry* topGrid, TopGridData* metaData);
	grid* first();
	grid* next();
};

struct SiblingIterator
{
	const HierarchyEntry* const firstSibling;
	const HierarchyEntry* currentSibling;

	SiblingIterator(const HierarchyEntry* const firstSibling);
	SiblingIterator(const LevelHierarchyEntry* const firstSibling);
	static SiblingIterator NewFromParent(const HierarchyEntry* const parent);
	static SiblingIterator NewFromParent(const LevelHierarchyEntry* const parent);

	grid* first();
	grid* next();
};

#endif /* SRC_ENZO_LEVELARRAYITERATOR_H_ */
