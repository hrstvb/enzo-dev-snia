cdef extern from "../enzo/Grid.h":
    # First we declare our class as being exposed to Cython
    ctypedef struct c_grid "grid":
        void DeleteAllFields()
        void AllocateGrids()
        int AppendForcingToBaryonFields() 
        int DetachForcingFromBaryonFields() 
        int SetExternalBoundaryValues(c_ExternalBoundary *Exterior)
    c_grid *new_grid "new grid" ()
    void del_grid "delete" (c_grid *g)

# Now we expose it to Python

cdef class grid:
    cdef c_grid *thisptr
    cdef int own
    def __cinit__(self, int own):
        self.thisptr = NULL
        self.own = own
        if self.own: self.thisptr = new_grid()
    def __dealloc__(self):
        if self.own: del_grid(self.thisptr)
    def DeleteAllFields(self):
        self.thisptr.DeleteAllFields()
    def AllocateGrids(self):
        self.thisptr.AllocateGrids()
    def AppendForcingToBaryonFields(self):
        self.thisptr.AppendForcingToBaryonFields()
    def DetachForcingFromBaryonFields(self):
        self.thisptr.DetachForcingFromBaryonFields()
    def SetExternalBoundaryValues(self, ExternalBoundary Exterior):
        self.thisptr.SetExternalBoundaryValues(Exterior.thisptr)
    def IsNull(self):
        return self.thisptr == NULL
