/// @file      newgrav-performance.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     lcaperf-related defines

#ifndef AMRSOLVE_PERFORMANCE_H
#define AMRSOLVE_PERFORMANCE_H

#ifdef USE_LCAPERF

#   include "lcaperf.h"

#   define LCAPERF_BEGIN(SEGMENT) lcaperf.begin(SEGMENT)
#   define LCAPERF_END(SEGMENT)   lcaperf.end(SEGMENT)
#   define LCAPERF_START(REGION)  lcaperf.start(REGION)
#   define LCAPERF_STOP(REGION)   lcaperf.stop(REGION)
#   define LCAPERF_GLOBAL(NAME,VALUE) lcaperf.global(NAME,VALUE)

#else /* USE_LCAPERF */

#   define LCAPERF_BEGIN(SEGMENT)     ;
#   define LCAPERF_END(SEGMENT)       ;
#   define LCAPERF_START(REGION)      ;
#   define LCAPERF_STOP(REGION)       ;
#   define LCAPERF_GLOBAL(NAME,VALUE) ;

#endif /* USE_LCAPERF */

#endif /* AMRSOLVE_PERFORMANCE_H */
