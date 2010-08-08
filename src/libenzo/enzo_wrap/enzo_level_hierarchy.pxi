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


