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

ctypedef float Eflt
ctypedef int Eint
ctypedef double FLOAT

cdef extern from "../enzo/typedefs.h":
    pass

include "enzo_globals.pxi"
global_data = _global_data()

cdef extern from "../enzo/Fluxes.h":
    pass

cdef extern from "../enzo/GridList.h":
    pass

cdef extern from "../enzo/ExternalBoundary.h":
    ctypedef struct c_ExternalBoundary "ExternalBoundary":
        pass
    c_ExternalBoundary *new_ExternalBoundary "new ExternalBoundary" ()
    void del_ExternalBoundary "delete" (c_ExternalBoundary *eb)

cdef class ExternalBoundary:
    cdef c_ExternalBoundary *thisptr
    def __cinit__(self):
        self.thisptr = new_ExternalBoundary()
    def __dealloc__(self):
        del_ExternalBoundary(self.thisptr)

cdef extern from "../enzo/Grid.h":
    # First we declare our class as being exposed to Cython
    ctypedef struct c_grid "grid":
        void DeleteAllFields()
        void AllocateGrids()
    c_grid *new_grid "new grid" ()
    void del_grid "delete" (c_grid *g)

# Now we expose it to Python

cdef class grid:
    cdef c_grid *thisptr
    def __cinit__(self):
        self.thisptr = new_grid()
    def __dealloc__(self):
        del_grid(self.thisptr)
    def DeleteAllFields(self):
        self.thisptr.DeleteAllFields()
    def AllocateGrids(self):
        self.thisptr.AllocateGrids()

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

cdef extern from "../enzo/TopGridData.h":
    struct c_TopGridData "TopGridData":
        pass

cdef class TopGridData:
    cdef c_TopGridData thisptr
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass

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

cdef class LevelHierarchyArray:
    cdef c_LevelHierarchyEntry *thisarray[MAX_DEPTH_OF_HIERARCHY]
    def __cinit__(self):
        cdef int i
        for i in range(MAX_DEPTH_OF_HIERARCHY):
            self.thisarray[i] = NULL
    def __dealloc__(self):
        pass
    def __getitem__(self, int i):
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
    int c_SetDefaultGlobalValues "SetDefaultGlobalValues" (c_TopGridData MetaData)
    int Group_ReadAllData(char *filename, c_HierarchyEntry *TopGrid, c_TopGridData tgd,
                          c_ExternalBoundary *Exterior)
    int c_CommunicationInitialize "CommunicationInitialize" (Eint32 *argc, char **argv[])
    int c_CommunicationPartitionGrid "CommunicationPartitionGrid" (c_HierarchyEntry *Grid, Eint gridnum)
    void c_AddLevel "AddLevel" (c_LevelHierarchyEntry *Array[], c_HierarchyEntry *Grid, Eint level)
    int c_EvolveHierarchy "EvolveHierarchy" (c_HierarchyEntry TopGrid, c_TopGridData tgd,
                        c_ExternalBoundary *Exterior, c_LevelHierarchyEntry *Array[],
                        Eflt Initialdt)

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

def read_all_data(char *filename, HierarchyEntry TopGrid, TopGridData MetaData, 
                  ExternalBoundary Exterior):
    Group_ReadAllData(filename, &TopGrid.thisptr, MetaData.thisptr,
                      Exterior.thisptr)

def AddLevel(LevelHierarchyArray Array, HierarchyEntry Grid, Eint level):
    c_AddLevel(Array.thisarray, &Grid.thisptr, level)

def EvolveHierarchy(HierarchyEntry TopGrid, TopGridData tgd,
                    ExternalBoundary Exterior, LevelHierarchyArray Array,
                    Eflt Initialdt):
    c_EvolveHierarchy(TopGrid.thisptr, tgd.thisptr, Exterior.thisptr,
                      Array.thisarray, Initialdt)

# This has to go at the end; it fixes all the "#define float double" stuff.
cdef extern from "fix_enzo_defs.h":
    pass
