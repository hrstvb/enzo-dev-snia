from stdlib cimport *

cdef extern from "math.h":
    pass

cdef extern from "mpi.h":
    pass

cdef extern from "../enzo/performance.h":
    pass

#cdef extern from "../enzo/ErrorExceptions.h":
    #pass

include "enzo_magic_numbers.pxi"

ctypedef double Eflt
ctypedef int Eint
ctypedef long int long_int
ctypedef double FLOAT

cdef extern from "../enzo/typedefs.h":
    pass

include "enzo_globals.pxi"
global_data = _global_data()

include "enzo_fluxes.pxi"

cdef extern from "../enzo/GridList.h":
    pass

include "enzo_external_boundary.pxi"
include "enzo_grid.pxi"

cdef extern from "../enzo/Hierarchy.h":
    cdef struct c_HierarchyEntry "HierarchyEntry":
        c_HierarchyEntry *NextGridThisLevel
        c_HierarchyEntry *NextGridNextLevel
        c_HierarchyEntry *ParentGrid
        c_grid         *GridData

cdef class HierarchyEntry:
    cdef c_HierarchyEntry thisptr
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass

include "enzo_top_grid_data.pxi"

cdef extern from "../enzo/LevelHierarchy.h":
    cdef struct c_LevelHierarchyEntry "LevelHierarchyEntry":
        c_LevelHierarchyEntry *NextGridThisLevel
        c_grid         *GridData
        c_HierarchyEntry *GridHierarchyEntry

cdef class LevelHierarchyEntry:
    cdef c_LevelHierarchyEntry *thisptr
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass

    property NextGridThisLevel:
        def __get__(self):
            cdef c_LevelHierarchyEntry *NGTL
            NGTL = self.thisptr.NextGridThisLevel
            if NGTL == NULL:
                return None
            cdef LevelHierarchyEntry lh = LevelHierarchyEntry()
            lh.thisptr = <c_LevelHierarchyEntry *> NGTL
            return lh

    property GridData:
        def __get__(self):
            cdef c_grid *GridData = self.thisptr.GridData
            cdef grid g = grid(False)
            g.thisptr = GridData
            return g

cdef class LevelHierarchyArray:
    cdef c_LevelHierarchyEntry *thisarray[MAX_DEPTH_OF_HIERARCHY]
    def __cinit__(self):
        cdef int i
        for i in range(MAX_DEPTH_OF_HIERARCHY):
            self.thisarray[i] = NULL
    def __dealloc__(self):
        pass
    def __getitem__(self, int i):
        #if self.thisarray[i] == NULL: return None
        cdef LevelHierarchyEntry a = LevelHierarchyEntry()
        a.thisptr = self.thisarray[i]
        return a

cdef extern from "../enzo/CommunicationUtilities.h":
    pass
 
cdef extern from "../enzo/macros_and_parameters.h":
    pass

ctypedef int Eint32

cdef extern from "../enzo/function_declarations.h":
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
    return c_SetDefaultGlobalValues(MetaData.thisptr)

def CommunicationPartitionGrid(HierarchyEntry Grid, int gridnum):
    return c_CommunicationPartitionGrid(&Grid.thisptr, gridnum)

def InitializeNew(char *filename, HierarchyEntry TopGrid, TopGridData MetaData,
                  ExternalBoundary Exterior):
    cdef Eflt Initialdt = 0.0
    c_InitializeNew(filename, TopGrid.thisptr, MetaData.thisptr,
                      Exterior.thisptr[0], &Initialdt)
    return Initialdt

def Group_ReadAllData(char *filename, HierarchyEntry TopGrid, TopGridData MetaData, 
                  ExternalBoundary Exterior):
    c_Group_ReadAllData(filename, &TopGrid.thisptr, MetaData.thisptr,
                      Exterior.thisptr)

def AddLevel(LevelHierarchyArray Array, HierarchyEntry Grid, Eint level):
    c_AddLevel(Array.thisarray, &Grid.thisptr, level)

def EvolveHierarchy(HierarchyEntry TopGrid, TopGridData tgd,
                    ExternalBoundary Exterior, LevelHierarchyArray Array,
                    Eflt Initialdt):
    c_EvolveHierarchy(TopGrid.thisptr, tgd.thisptr, Exterior.thisptr,
                      Array.thisarray, Initialdt)

def CopyOverlappingZones(grid CurrentGrid, TopGridData MetaData,
                         LevelHierarchyArray Array, int level):
    cdef int rv
    rv = c_CopyOverlappingZones(<c_grid *> grid.thisptr, &MetaData.thisptr,
                                Array.thisarray, level)
    return rv

def CommunicationReceiveHandler():
    c_CommunicationReceiveHandler()

# This has to go at the end; it fixes all the "#define float double" stuff.
cdef extern from "fix_enzo_defs.h":
    pass
