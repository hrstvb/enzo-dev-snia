// $Id: newgrav-performance.h 10 2010-03-18 22:35:19Z jbordner $

/// @file      newgrav-performance.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     lcaperf-related defines

#ifndef NEWGRAV_PERFORMANC_H
#define NEWGRAV_PERFORMANC_H

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

#endif /* NEWGRAV_PERFORMANC_H */
