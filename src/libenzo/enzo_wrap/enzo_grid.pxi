cdef extern from "Grid.h":
    # First we declare our class as being exposed to Cython
    ctypedef struct c_grid "grid":
        void DeleteAllFields()
        void AllocateGrids()
        Eint AppendForcingToBaryonFields() 
        Eint DetachForcingFromBaryonFields() 
        Eint SetExternalBoundaryValues(c_ExternalBoundary *Exterior)

        # All the grid properties go here
        
        Eint GridRank
        Eint GridDimension[MAX_DIMENSION]
        Eint GridStartIndex[MAX_DIMENSION]
        Eint GridEndIndex[MAX_DIMENSION]

        Eflt   GridLeftEdge[MAX_DIMENSION]
        Eflt   GridRightEdge[MAX_DIMENSION]
        Eflt   dtFixed
        Eflt   Time
        Eflt   OldTime
        Eint   SubgridsAreStatic
        Eint   ID
        
        # Baryon grid data
        
        Eint    NumberOfBaryonFields
        Eflt   *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS]
        Eflt   *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]
        Eflt   *InterpolatedField[MAX_NUMBER_OF_BARYON_FIELDS]
        Eflt   *RandomForcingField[MAX_DIMENSION]
        Eint    FieldType[MAX_NUMBER_OF_BARYON_FIELDS]
        Eflt   *CellLeftEdge[MAX_DIMENSION]
        Eflt   *CellWidth[MAX_DIMENSION]
        c_fluxes *BoundaryFluxes

        # MHD data

        Eflt   *divB
        Eflt   *gradPhi[MAX_DIMENSION]

        Eflt    CourantSafetyNumber
        Eint    PPMFlatteningParameter
        Eint    PPMDiffusionParameter
        Eint    PPMSteepeningParameter
        
        # Particle data
        
        Eint    NumberOfParticles
        Eflt   *ParticlePosition[MAX_DIMENSION]
        Eflt   *ParticleVelocity[MAX_DIMENSION]
        Eflt   *ParticleAcceleration[MAX_DIMENSION+1]
        Eflt   *ParticleMass
        Eint   *ParticleNumber
        Eint   *ParticleType
        Eflt   *ParticleAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES]
        
        # Star particle data
        
        Eint NumberOfStars
        #Star *Stars

        # Gravity data
     
        Eflt   *PotentialField
        Eflt   *AccelerationField[MAX_DIMENSION]
        Eflt   *GravitatingMassField
        Eflt    GravitatingMassFieldLeftEdge[MAX_DIMENSION]
        Eint    GravitatingMassFieldDimension[MAX_DIMENSION]
        Eflt    GravitatingMassFieldCellSize
        Eflt   *GravitatingMassFieldParticles
        Eflt    GravitatingMassFieldParticlesLeftEdge[MAX_DIMENSION]
        Eflt    GravitatingMassFieldParticlesCellSize
        Eint    GravitatingMassFieldParticlesDimension[MAX_DIMENSION]

        #gravity_boundary_type GravityBoundaryType

        Eflt    PotentialSum
        
        # Rebuild Hierarchy Temporaries
        
        Eint   *FlaggingField
        Eflt   *MassFlaggingField
        Eflt   *ParticleMassFlaggingField
        
        # Parallel Information
        
        Eint ProcessorNumber

        #  Movie Data Format
        Eint TimestepsSinceCreation

    c_grid *new_grid "new grid" ()
    void del_grid "delete" (c_grid *g)

# Now we expose it to Python

cdef class grid:
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

    def get_baryon_field(self, int index, int old = 0):
        cdef int dim, size = 1
        cdef Eint32 ndims = self.thisptr.GridRank
        cdef np.npy_intp dims[3]
        for dim in range(ndims):
            size *= self.thisptr.GridDimension[dim]
            dims[dim] = self.thisptr.GridDimension[dim]
        cdef Eflt *bfield
        if old: bfield = self.thisptr.OldBaryonField[index]
        else: bfield = self.thisptr.BaryonField[index]
        if bfield == NULL: return None
        cdef np.ndarray field
        field = cnp.PyArray_SimpleNewFromData(
                    ndims, dims, E_ENPY_BFLOAT,
                    bfield)
        # Might need to set field flags here
        return field

    def get_flux(self, int index):
        cdef c_grid *g = self.thisptr
        fluxes = []
        # We are essentially replicating Grid_ClearBoundaryFluxes's math
        cdef Eint32 ndims = self.thisptr.GridRank
        cdef int dim, i, field
        cdef np.npy_intp size
        cdef np.npy_intp fdims = 1
        cdef np.ndarray f1, f2
        for dim in range(ndims):
            size = 1
            for i in range(ndims):
                size *= (g.BoundaryFluxes.LeftFluxEndGlobalIndex[dim][i] -
                         g.BoundaryFluxes.LeftFluxStartGlobalIndex[dim][i] + 1)
            for field in range(g.NumberOfBaryonFields):
                if g.BoundaryFluxes.LeftFluxes[field][dim] == NULL:
                    f1 = None
                else:
                    f1 = cnp.PyArray_SimpleNewFromData(
                                fdims, &size, E_ENPY_BFLOAT,
                                g.BoundaryFluxes.LeftFluxes[field][dim])
                if g.BoundaryFluxes.RightFluxes[field][dim] == NULL:
                    f2 = None
                else:
                    f2 = cnp.PyArray_SimpleNewFromData(
                                fdims, &size, E_ENPY_BFLOAT,
                                g.BoundaryFluxes.RightFluxes[field][dim])
                fluxes.append((f1,f2))
        return fluxes

    property ParticlePositions:
        def __get__(self):
            cdef int dim
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef FLOAT **PField = self.thisptr.ParticlePosition
            pos = []
            for dim in range(MAX_DIMENSION):
                if PField[dim] == NULL:
                    pos.append(None)
                    continue
                pos.append(cnp.PyArray_SimpleNewFromData(
                               1, &NumberOfParticles, E_ENPY_PFLOAT,
                               PField[dim]))
            return pos

    property ParticleVelocity:
        def __get__(self):
            cdef int dim
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef Eflt **PField = self.thisptr.ParticleVelocity
            pos = []
            for dim in range(MAX_DIMENSION):
                if PField[dim] == NULL:
                    pos.append(None)
                    continue
                pos.append(cnp.PyArray_SimpleNewFromData(
                               1, &NumberOfParticles, E_ENPY_BFLOAT,
                               PField[dim]))
            return pos

    property ParticleAcceleration:
        def __get__(self):
            cdef int dim
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef Eflt **PField = self.thisptr.ParticleAcceleration
            pos = []
            for dim in range(MAX_DIMENSION+1):
                if PField[dim] == NULL:
                    pos.append(None)
                    continue
                pos.append(cnp.PyArray_SimpleNewFromData(
                               1, &NumberOfParticles, E_ENPY_BFLOAT,
                               PField[dim]))
            return pos

    property ParticleMass:
        def __get__(self):
            cdef np.ndarray f
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef Eflt *PField = self.thisptr.ParticleMass
            if PField == NULL:
                return None
            f = cnp.PyArray_SimpleNewFromData(
                           1, &NumberOfParticles, E_ENPY_BFLOAT,
                           PField)
            return f

    property ParticleNumber:
        def __get__(self):
            cdef np.ndarray f
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef Eint *PField = self.thisptr.ParticleNumber
            if PField == NULL:
                return None
            f = cnp.PyArray_SimpleNewFromData(
                           1, &NumberOfParticles, E_ENPY_INT,
                           PField)
            return f

    property ParticleType:
        def __get__(self):
            cdef np.ndarray f
            cdef np.npy_intp TypeOfParticles = self.thisptr.NumberOfParticles
            cdef Eint *PField = self.thisptr.ParticleType
            if PField == NULL:
                return None
            f = cnp.PyArray_SimpleNewFromData(
                           1, &TypeOfParticles, E_ENPY_INT,
                           PField)
            return f

    property ParticleAttribute:
        def __get__(self):
            cdef int dim
            cdef np.npy_intp NumberOfParticles = self.thisptr.NumberOfParticles
            cdef Eflt **PField = self.thisptr.ParticleAttribute
            pos = []
            for dim in range(MAX_NUMBER_OF_PARTICLE_ATTRIBUTES):
                if PField[dim] == NULL:
                    pos.append(None)
                    continue
                pos.append(cnp.PyArray_SimpleNewFromData(
                               1, &NumberOfParticles, E_ENPY_BFLOAT,
                               PField[dim]))
            return pos

    # This code was auto-generated by wrap_dot_h.py

    property GridRank:
        def __get__(self):
            return self.thisptr.GridRank
        def __set__(self, Eint val):
            self.thisptr.GridRank = val


    property GridDimension:
        def __get__(self):
            print "Returning a copy of GridDimension"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GridDimension[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GridDimension[i] = val[i]


    property GridStartIndex:
        def __get__(self):
            print "Returning a copy of GridStartIndex"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GridStartIndex[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GridStartIndex[i] = val[i]


    property GridEndIndex:
        def __get__(self):
            print "Returning a copy of GridEndIndex"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GridEndIndex[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GridEndIndex[i] = val[i]


    property GridLeftEdge:
        def __get__(self):
            print "Returning a copy of GridLeftEdge"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GridLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GridLeftEdge[i] = val[i]


    property GridRightEdge:
        def __get__(self):
            print "Returning a copy of GridRightEdge"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GridRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GridRightEdge[i] = val[i]


    property dtFixed:
        def __get__(self):
            return self.thisptr.dtFixed
        def __set__(self, Eflt val):
            self.thisptr.dtFixed = val


    property Time:
        def __get__(self):
            return self.thisptr.Time
        def __set__(self, Eflt val):
            self.thisptr.Time = val


    property OldTime:
        def __get__(self):
            return self.thisptr.OldTime
        def __set__(self, Eflt val):
            self.thisptr.OldTime = val


    property SubgridsAreStatic:
        def __get__(self):
            return self.thisptr.SubgridsAreStatic
        def __set__(self, Eint val):
            self.thisptr.SubgridsAreStatic = val


    property ID:
        def __get__(self):
            return self.thisptr.ID
        def __set__(self, Eint val):
            self.thisptr.ID = val


    property NumberOfBaryonFields:
        def __get__(self):
            return self.thisptr.NumberOfBaryonFields
        def __set__(self, Eint val):
            self.thisptr.NumberOfBaryonFields = val


    property FieldType:
        def __get__(self):
            print "Returning a copy of FieldType"
            retval = []
            for i in range(MAX_NUMBER_OF_BARYON_FIELDS):
                retval.append(self.thisptr.FieldType[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_NUMBER_OF_BARYON_FIELDS):
                self.thisptr.FieldType[i] = val[i]


    property CourantSafetyNumber:
        def __get__(self):
            return self.thisptr.CourantSafetyNumber
        def __set__(self, Eflt val):
            self.thisptr.CourantSafetyNumber = val


    property PPMFlatteningParameter:
        def __get__(self):
            return self.thisptr.PPMFlatteningParameter
        def __set__(self, Eint val):
            self.thisptr.PPMFlatteningParameter = val


    property PPMDiffusionParameter:
        def __get__(self):
            return self.thisptr.PPMDiffusionParameter
        def __set__(self, Eint val):
            self.thisptr.PPMDiffusionParameter = val


    property PPMSteepeningParameter:
        def __get__(self):
            return self.thisptr.PPMSteepeningParameter
        def __set__(self, Eint val):
            self.thisptr.PPMSteepeningParameter = val


    property NumberOfParticles:
        def __get__(self):
            return self.thisptr.NumberOfParticles
        def __set__(self, Eint val):
            self.thisptr.NumberOfParticles = val


    property NumberOfStars:
        def __get__(self):
            return self.thisptr.NumberOfStars
        def __set__(self, Eint val):
            self.thisptr.NumberOfStars = val


    property GravitatingMassFieldLeftEdge:
        def __get__(self):
            print "Returning a copy of GravitatingMassFieldLeftEdge"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GravitatingMassFieldLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GravitatingMassFieldLeftEdge[i] = val[i]


    property GravitatingMassFieldDimension:
        def __get__(self):
            print "Returning a copy of GravitatingMassFieldDimension"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GravitatingMassFieldDimension[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GravitatingMassFieldDimension[i] = val[i]


    property GravitatingMassFieldCellSize:
        def __get__(self):
            return self.thisptr.GravitatingMassFieldCellSize
        def __set__(self, Eflt val):
            self.thisptr.GravitatingMassFieldCellSize = val


    property GravitatingMassFieldParticlesLeftEdge:
        def __get__(self):
            print "Returning a copy of GravitatingMassFieldParticlesLeftEdge"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GravitatingMassFieldParticlesLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GravitatingMassFieldParticlesLeftEdge[i] = val[i]


    property GravitatingMassFieldParticlesCellSize:
        def __get__(self):
            return self.thisptr.GravitatingMassFieldParticlesCellSize
        def __set__(self, Eflt val):
            self.thisptr.GravitatingMassFieldParticlesCellSize = val


    property GravitatingMassFieldParticlesDimension:
        def __get__(self):
            print "Returning a copy of GravitatingMassFieldParticlesDimension"
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(self.thisptr.GravitatingMassFieldParticlesDimension[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(MAX_DIMENSION):
                self.thisptr.GravitatingMassFieldParticlesDimension[i] = val[i]


    property PotentialSum:
        def __get__(self):
            return self.thisptr.PotentialSum
        def __set__(self, Eflt val):
            self.thisptr.PotentialSum = val


    property ProcessorNumber:
        def __get__(self):
            return self.thisptr.ProcessorNumber
        def __set__(self, Eint val):
            self.thisptr.ProcessorNumber = val


    property TimestepsSinceCreation:
        def __get__(self):
            return self.thisptr.TimestepsSinceCreation
        def __set__(self, Eint val):
            self.thisptr.TimestepsSinceCreation = val




# SKIPPED and dealt with:
# SKIPPED: '*BaryonField[MAX_NUMBER_OF_BARYON_FIELDS]'
# SKIPPED: '*OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]'

# SKIPPED and not dealt with:

# SKIPPED: '*InterpolatedField[MAX_NUMBER_OF_BARYON_FIELDS]'
# SKIPPED: '*RandomForcingField[MAX_DIMENSION]'
# SKIPPED: '*CellLeftEdge[MAX_DIMENSION]'
# SKIPPED: '*CellWidth[MAX_DIMENSION]'
# SKIPPED: '*divB'
# SKIPPED: '*gradPhi[MAX_DIMENSION]'
# SKIPPED: '*ParticlePosition[MAX_DIMENSION]'
# SKIPPED: '*ParticleVelocity[MAX_DIMENSION]'
# SKIPPED: '*ParticleAcceleration[MAX_DIMENSION+1]'
# SKIPPED: '*ParticleMass'
# SKIPPED: '*ParticleNumber'
# SKIPPED: '*ParticleType'
# SKIPPED: '*ParticleAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES]'
# SKIPPED: '*PotentialField'
# SKIPPED: '*AccelerationField[MAX_DIMENSION]'
# SKIPPED: '*GravitatingMassField'
# SKIPPED: '*GravitatingMassFieldParticles'
# SKIPPED: '*FlaggingField'
# SKIPPED: '*MassFlaggingField'
# SKIPPED: '*ParticleMassFlaggingField'
