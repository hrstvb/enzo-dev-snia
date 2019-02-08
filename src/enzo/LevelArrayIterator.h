/*
 * LevelArrayIterator.h
 *
 *  Created on: Jan 23, 2019
 *      Author: cdev
 */

#ifndef SRC_ENZO_LEVELARRAYITERATOR_H_
#define SRC_ENZO_LEVELARRAYITERATOR_H_

class grid;
#include "LevelHierarchy.h"

struct LevelArrayIterator
{
	LevelHierarchyEntry** levelArray; //
	int nLevels = -1; ///
	bool forward = true; //
	int currentLevel = 0; //
	LevelHierarchyEntry* currentEntry = NULL;

	LevelArrayIterator(LevelHierarchyEntry** levelArray);
	LevelArrayIterator(LevelHierarchyEntry** levelArray, int numberOfLevels);int countLevels();
	void resetToTop();
	void resetToFinest();
	grid* firstFromTop();
	grid* firstFromTop(grid** parent);
	grid* firstFromTop(int* level);
	grid* firstFromTop(int* level, grid** parent);
	grid* firstFromFinest();
	grid* firstFromFinest(grid** parent);
	grid* firstFromFinest(int* level);
	grid* firstFromFinest(int* level, grid** parent);
	grid* next();
	grid* next(grid** parent);
	grid* next(int* level);
	grid* next(int* level, grid** parent);
};

#endif /* SRC_ENZO_LEVELARRAYITERATOR_H_ */
