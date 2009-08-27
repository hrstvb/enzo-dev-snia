ctypedef double Eflt
ctypedef int Eint
ctypedef long int long_int
ctypedef double FLOAT
ctypedef int Eint32
ctypedef int Eint64

# All our incldues...

cdef extern from "math.h":
    pass

cdef extern from "mpi.h":
    pass

cdef extern from "performance.h":
    pass

cdef extern from "macros_and_parameters.h":
    pass

cdef extern from "typedefs.h":
    pass

cdef extern from "global_data.h":
    pass

cdef extern from "Fluxes.h":
    pass

cdef extern from "GridList.h":
    pass

cdef extern from "ExternalBoundary.h":
    ctypedef struct c_ExternalBoundary "ExternalBoundary"

cdef extern from "Grid.h":
    ctypedef struct c_grid "grid"

cdef extern from "Hierarchy.h":
    cdef struct c_HierarchyEntry "HierarchyEntry":
        c_HierarchyEntry *NextGridThisLevel
        c_HierarchyEntry *NextGridNextLevel
        c_HierarchyEntry *ParentGrid
        c_grid         *GridData

cdef extern from "TopGridData.h":
    cdef struct c_TopGridData "TopGridData"

cdef extern from "LevelHierarchy.h":
    cdef struct c_LevelHierarchyEntry "LevelHierarchyEntry":
        c_LevelHierarchyEntry *NextGridThisLevel
        c_grid         *GridData
        c_HierarchyEntry *GridHierarchyEntry

cdef extern from "CommunicationUtilities.h":
    pass

# Now we forward declare all our cdef classes
 
cdef class ExternalBoundary:
    cdef c_ExternalBoundary *thisptr

cdef class grid:
    cdef c_grid *thisptr
    cdef int own

cdef class TopGridData:
    cdef c_TopGridData *thisptr
