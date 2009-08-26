cdef class grid

cdef extern from "../enzo/ExternalBoundary.h":
    ctypedef struct c_grid "grid::grid" # We have to namespace this!
    ctypedef struct c_ExternalBoundary "ExternalBoundary":
        # Methods on the class
        int AppendForcingToBaryonFields()
        int DetachForcingFromBaryonFields()
        int Prepare(c_grid *TopGrid)

    # Allocation / deallocation
    c_ExternalBoundary *new_ExternalBoundary "new ExternalBoundary" ()
    void del_ExternalBoundary "delete" (c_ExternalBoundary *eb)

cdef class ExternalBoundary:
    cdef c_ExternalBoundary *thisptr
    def __cinit__(self):
        self.thisptr = new_ExternalBoundary()
    def __dealloc__(self):
        del_ExternalBoundary(self.thisptr)
    def AppendForcingToBaryonFields(self):
        self.thisptr.AppendForcingToBaryonFields()
    def DetachForcingFromBaryonFields(self):
        self.thisptr.DetachForcingFromBaryonFields()
    def Prepare(self, grid TopGrid):
        cdef c_grid *cg = <c_grid *> grid.thisptr
        self.thisptr.Prepare(cg)
