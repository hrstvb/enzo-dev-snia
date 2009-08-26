Eint32 enzo_main(Eint32 argc, char** argv);
int SetDefaultGlobalValues(TopGridData &MetaData);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    LevelHierarchyEntry *Array[], float Initialdt);

