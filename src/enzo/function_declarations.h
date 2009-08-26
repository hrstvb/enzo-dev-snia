Eint32 enzo_main(Eint32 argc, char** argv);
int SetDefaultGlobalValues(TopGridData &MetaData);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
