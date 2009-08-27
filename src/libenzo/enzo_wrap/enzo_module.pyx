from stdlib cimport *

include "enzo_magic_numbers.pxi"
include "enzo_globals.pxi"
include "enzo_fluxes.pxi"
include "enzo_external_boundary.pxi"
include "enzo_grid.pxi"
include "enzo_hierarchy.pxi"
include "enzo_top_grid_data.pxi"
include "enzo_level_hierarchy.pxi"

global_data = _global_data()

cdef extern from "function_declarations.h":
    # All functions that don't get prototyped elsewhere get prototyped here.
    # This fixes name mangling issues.
    Eint32 c_enzo_main "enzo_main" (Eint32 argc, char **argv) except +
    int c_SetDefaultGlobalValues "SetDefaultGlobalValues" (c_TopGridData MetaData) except +
    int c_InitializeNew "InitializeNew" (char *filename, c_HierarchyEntry TopGrid,
            c_TopGridData tgd, c_ExternalBoundary Exterior, Eflt *Initialdt) except +
    int c_Group_ReadAllData "Group_ReadAllData" (
                          char *filename, c_HierarchyEntry *TopGrid, c_TopGridData tgd,
                          c_ExternalBoundary *Exterior) except +
    int c_CommunicationInitialize "CommunicationInitialize" (Eint32 *argc, char **argv[]) except +
    int c_CommunicationPartitionGrid "CommunicationPartitionGrid" (c_HierarchyEntry *Grid, Eint gridnum)
    void c_AddLevel "AddLevel" (c_LevelHierarchyEntry *Array[], c_HierarchyEntry *Grid, Eint level)
    int c_EvolveHierarchy "EvolveHierarchy" (c_HierarchyEntry TopGrid, c_TopGridData tgd,
                        c_ExternalBoundary *Exterior, c_LevelHierarchyEntry *Array[],
                        Eflt Initialdt) except +
    int c_CopyOverlappingZones "CopyOverlappingZones" (
                             c_grid* CurrentGrid, c_TopGridData *MetaData,
                             c_LevelHierarchyEntry *LevelArray[], int level)
    int c_CommunicationReceiveHandler "CommunicationReceiveHandler" (
            c_fluxes **SubGridFluxesEstimate[] = NULL,
            int NumberOfSubgrids[] = NULL, int FluxFlag = FALSE,
            c_TopGridData *MetaData = NULL)
    int c_CommunicationBarrier "CommunicationBarrier" ()
    int c_EvolveLevel "EvolveLevel" (c_TopGridData *MetaData, c_LevelHierarchyEntry *Array[],
                        int level, Eflt dtLevelAbove, c_ExternalBoundary *Exterior) except +

def run_enzo_main(args):
    cdef int argc = len(args)
    cdef char **argv = <char **>malloc(argc * sizeof(char*))
    for i in range(argc):
        argv[i] = <char *>args[i]
    c_enzo_main(argc, argv)
    free(argv)

def CommunicationInitialize(args):
    cdef int argc = len(args)
    cdef int *argcs = &argc
    cdef char **argv = <char **>malloc(argc * sizeof(char*))
    cdef char ***argvs = &argv
    for i in range(argc):
        argv[i] = <char *>args[i]
    c_CommunicationInitialize(argcs, argvs)
    free(argv)

def SetDefaultGlobalValues(TopGridData MetaData):
    return c_SetDefaultGlobalValues(MetaData.thisptr[0])

def CommunicationPartitionGrid(HierarchyEntry Grid, int gridnum):
    return c_CommunicationPartitionGrid(&Grid.thisptr, gridnum)

def InitializeNew(char *filename, HierarchyEntry TopGrid, TopGridData MetaData,
                  ExternalBoundary Exterior):
    cdef Eflt Initialdt = 0.0
    c_InitializeNew(filename, TopGrid.thisptr, MetaData.thisptr[0],
                      Exterior.thisptr[0], &Initialdt)
    return Initialdt

def Group_ReadAllData(char *filename, HierarchyEntry TopGrid, TopGridData MetaData, 
                  ExternalBoundary Exterior):
    c_Group_ReadAllData(filename, &TopGrid.thisptr, MetaData.thisptr[0],
                      Exterior.thisptr)

def AddLevel(LevelHierarchyArray Array, HierarchyEntry Grid, Eint level):
    c_AddLevel(Array.thisarray, &Grid.thisptr, level)

def EvolveHierarchy(HierarchyEntry TopGrid, TopGridData MetaData,
                    ExternalBoundary Exterior, LevelHierarchyArray Array,
                    Eflt Initialdt):
    c_EvolveHierarchy(TopGrid.thisptr, MetaData.thisptr[0], Exterior.thisptr,
                      Array.thisarray, Initialdt)

def CopyOverlappingZones(grid CurrentGrid, TopGridData MetaData,
                         LevelHierarchyArray Array, int level):
    cdef int rv
    rv = c_CopyOverlappingZones(<c_grid *> CurrentGrid.thisptr, MetaData.thisptr,
                                Array.thisarray, level)
    return rv

def CommunicationReceiveHandler():
    c_CommunicationReceiveHandler()

def EvolveLevel(TopGridData MetaData, LevelHierarchyArray Array,
                int level, Eflt dtLevelAbove, ExternalBoundary Exterior):
    c_EvolveLevel(MetaData.thisptr, Array.thisarray,
                  level, dtLevelAbove, Exterior.thisptr)

# This has to go at the end; it fixes all the "#define float double" stuff.
cdef extern from "fix_enzo_defs.h":
    pass
