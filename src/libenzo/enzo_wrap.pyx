from stdlib cimport *

cdef extern from "math.h":
    pass

cdef extern from "mpi.h":
    pass

cdef extern from "../enzo/performance.h":
    pass

#cdef extern from "../enzo/ErrorExceptions.h":
    #pass

cdef extern from "../enzo/macros_and_parameters.h":
    pass

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
    cdef struct HierarchyEntry:
        HierarchyEntry *NextGridThisLevel
        HierarchyEntry *NextGridNextLevel
        HierarchyEntry *ParentGrid
        c_grid         *GridData

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
    pass

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
    int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, c_TopGridData tgd,
                          c_ExternalBoundary *Exterior)

def run_enzo_main(args):
    cdef int argc = len(args)
    cdef char **argv = <char **>malloc(argc * sizeof(char*))
    for i in range(argc):
        argv[i] = <char *>args[i]
    c_enzo_main(argc, argv)

def SetDefaultGlobalValues(TopGridData MetaData):
    return c_SetDefaultGlobalValues(MetaData.thisptr)

def read_all_data(char *filename, TopGridData MetaData):
    cdef HierarchyEntry *TopGrid = <HierarchyEntry *>malloc(sizeof(HierarchyEntry))
    cdef c_ExternalBoundary *Exterior = new_ExternalBoundary()
    Group_ReadAllData(filename, TopGrid, MetaData.thisptr, Exterior)

# This has to go at the end; it fixes all the "#define float double" stuff.
cdef extern from "fix_enzo_defs.h":
    pass
