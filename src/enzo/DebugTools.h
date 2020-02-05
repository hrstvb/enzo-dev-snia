#ifndef DEBUG_TOOLS_H
#define DEBUG_TOOLS_H

class grid;
struct TopGridData;
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

int sprintHierarchy(char* s, HierarchyEntry* topGrid, const char* const filename, const int linenum);
int printHierarchy(HierarchyEntry* topGrid);
int printHierarchy0(HierarchyEntry* topGrid);
int printHierarchy(HierarchyEntry* topGrid, const char* const filename, const int linenum);
int printHierarchy0(HierarchyEntry* topGrid, const char* const filename, const int linenum);
int sprintHierarchy(char* s, LevelHierarchyEntry** levelArray, const char* const filename, const int linenum);
int printHierarchy(LevelHierarchyEntry** levelArray);
int printHierarchy0(LevelHierarchyEntry** levelArray);
int printHierarchy(LevelHierarchyEntry** levelArray, const char* const filename, const int linenum);
int printHierarchy0(LevelHierarchyEntry** levelArray, const char* const filename, const int linenum);
#define PRINT_HIERARCHY(levelArray) printHierarchy(levelArray, __FILE__, __LINE__)
#define PRINT_HIERARCHY0(levelArray) printHierarchy0(levelArray, __FILE__, __LINE__)

#endif /* DEBUG_TOOLS_H */
