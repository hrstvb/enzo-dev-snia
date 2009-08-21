cdef extern from "../enzo/performance.h":
    pass

#cdef extern from "../enzo/ErrorExceptions.h":
    #pass

cdef extern from "../enzo/macros_and_parameters.h":
    pass

cdef extern from "../enzo/typedefs.h":
    pass

cdef extern from "../enzo/global_data.h":
    pass

cdef extern from "../enzo/Fluxes.h":
    pass

cdef extern from "../enzo/GridList.h":
    pass

cdef extern from "../enzo/ExternalBoundary.h":
    pass

cdef extern from "../enzo/Grid.h":
    ctypedef struct c_grid "grid":
        void DeleteAllFields()
        void AllocateGrids()
    c_grid *new_grid "new grid" ()
    void del_grid "delete grid" ()

cdef extern from "../enzo/Hierarchy.h":
    pass

cdef extern from "../enzo/TopGridData.h":
    pass

cdef extern from "../enzo/LevelHierarchy.h":
    pass

cdef extern from "../enzo/CommunicationUtilities.h":
    pass
 
cdef extern from "../enzo/macros_and_parameters.h":
    pass
