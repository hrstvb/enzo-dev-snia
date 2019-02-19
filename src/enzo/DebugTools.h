#ifndef DEBUG_TOOLS_H
#define DEBUG_TOOLS_H

class grid;
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "DebugMacros.h"

void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label);

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid,
				    TopGridData *MetaData);

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

size_t sprintHierarchy(char* s, LevelHierarchyEntry** levelArray);
void printHierarchy(LevelHierarchyEntry** levelArray);

#endif /* DEBUG_TOOLS_H */
