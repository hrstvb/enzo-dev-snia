#ifndef DEBUG_TOOLS_H_
#define DEBUG_TOOLS_H_

void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label);

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid,
				    TopGridData *MetaData);

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif


#define DEF2STR_(macro) #macro
#define DEF2STR(macro) DEF2STR_(macro)
#define DEF2FMT(macro) DEF2STR(__FILE__)":"DEF2STR(__LINE__)": "#macro "="  DEF2STR(macro)
#define PRAGMA_SHOW_DEF(macro) "message("DEF2FMT(macro)")"

#endif /* DEBUG_TOOLS_H_ */
