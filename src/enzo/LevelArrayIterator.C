#include "LevelArrayIterator.h"

LevelArrayIterator::LevelArrayIterator(LevelHierarchyEntry** levelArray)
{
	this->levelArray = levelArray;
}

LevelArrayIterator::LevelArrayIterator(LevelHierarchyEntry** levelArray, int numberOfLevels)
{
	this->levelArray = levelArray;
	nLevels = (numberOfLevels >= 0) ? numberOfLevels : countLevels();
}

int LevelArrayIterator::countLevels()
{
	int n = 0;
	for(; n <= MAX_DEPTH_OF_HIERARCHY; n++)
		if(levelArray[n] == NULL)
			break;
	return n;
}

void LevelArrayIterator::resetToTop()
{
	currentEntry = NULL;
	forward = true;
}

void LevelArrayIterator::resetToFinest()
{
	currentEntry = NULL;
	forward = false;
}

grid* LevelArrayIterator::firstFromTop()
{
	resetToTop();
	return next();
}

grid* LevelArrayIterator::firstFromTop(grid** parent)
{
	resetToTop();
	return next(parent);
}

grid* LevelArrayIterator::firstFromTop(int* level)
{
	resetToTop();
	return next(level);
}

grid* LevelArrayIterator::firstFromTop(int* level, grid** parent)
{
	resetToTop();
	return next(level, parent);
}

grid* LevelArrayIterator::firstFromFinest()
{
	resetToFinest();
	return next();
}

grid* LevelArrayIterator::firstFromFinest(grid** parent)
{
	resetToFinest();
	return next(parent);
}

grid* LevelArrayIterator::firstFromFinest(int* level)
{
	resetToFinest();
	return next(level);
}

grid* LevelArrayIterator::firstFromFinest(int* level, grid** parent)
{
	resetToFinest();
	return next(level, parent);
}

grid* LevelArrayIterator::next()
{
	int level;
	grid* parent;
	return next(&level, &parent);
}

grid* LevelArrayIterator::next(grid** parent)
{
	int level;
	return next(&level, parent);
}

grid* LevelArrayIterator::next(int* level)
{
	grid* parent;
	return next(level, &parent);
}

grid* LevelArrayIterator::next(int* level, grid** parent)
{
	if(currentEntry)
	{
		currentEntry = currentEntry->NextGridThisLevel;
		if(currentEntry == NULL)
			currentLevel += (forward) ? 1 : -1;
	}
	else
	{
		int nLevelsTemp = (nLevels >= 0) ? nLevels : countLevels();
		currentLevel = (nLevelsTemp > 0) ? ((forward) ? 0 : nLevelsTemp - 1) : -1;
	}

	*level = currentLevel;
	if(currentLevel < 0 || currentLevel > MAX_DEPTH_OF_HIERARCHY)
		return *parent = NULL;

	if(currentEntry == NULL)
	{
		currentEntry = levelArray[currentLevel];
		if(currentEntry == NULL)
			return *parent = NULL;
	}

	HierarchyEntry* he = (currentEntry) ? currentEntry->GridHierarchyEntry : NULL;
	he = (he) ? he->ParentGrid : he;
	*parent = (he) ? he->GridData : NULL;
	return currentEntry->GridData;
}
