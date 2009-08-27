cdef extern from "Fluxes.h":
    cdef struct c_fluxes "fluxes":
        long_int LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION]
        long_int LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION]
        long_int RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION]
        long_int RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION]
        Eflt *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION]
        Eflt *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION]
