#ifndef DEBUG_TOOLS_H
#define DEBUG_TOOLS_H

void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label);

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid,
				    TopGridData *MetaData);

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

#define MACRO2STRING_(x) #x
#define MACRO2STR(x) MACRO2STRING_(x)
#define MACRO_FMT(macro) MACRO2STR(__FILE__)":"MACRO2STR(__LINE__)": "#macro "="  MACRO2STR(macro)
#define PRAGMA_MACRO(macro) "message("MACRO_FMT(macro)")"

#endif /* DEBUG_TOOLS_H */
