#ifdef __macros_and_parameters_h_
ERROR: need performance.h to be included before macros_and_parameters.h
#endif

/*****************************************************************************
 *                                                                           *
 * Copyright 2006 James Bordner
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Board of Trustees of the University of Illinois            *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
#ifndef PERFORMANCE_H
#define PERFORMANCE_H

//======================================================================
//
// File:        performance.h
//
// Description: Interface layer between Enzo and lcaperf
//
// To use lcaperf, compile with -DUSE_LCAPERF and link with -llcaperf
// To use PAPI with lcaperf, compile also with -DUSE_PAPI and link with -lpapi
//
//----------------------------------------------------------------------
//
// James Bordner (jobordner@ucsd.edu)
// 2003-06-20
//
//======================================================================

//----------------------------------------------------------------------
// lcaperf
//----------------------------------------------------------------------
#ifdef USE_LCAPERF

#   include "lcaperf.hpp" // THIS MUST COME BEFORE int IS REDEFINED

#  define LCAPERF_BEGIN(segment)    lcaperf.begin ()
#  define LCAPERF_END(segment)      lcaperf.end ()
#  define LCAPERF_START(region)     lcaperf.start (region)
#  define LCAPERF_STOP(region)      lcaperf.stop (region)
#  define LCAPERF_PRINT             lcaperf.print()

class LcaPerfEnzo : public LcaPerf
{
 public:
  virtual void header ();
  virtual void print ();
};

extern LcaPerfEnzo lcaperf;

#else

#  define LCAPERF_BEGIN(segment)    /* This space intentionally left blank */ ;
#  define LCAPERF_END(segment)      /* This space intentionally left blank */ ;
#  define LCAPERF_START(region)     /* This space intentionally left blank */ ;
#  define LCAPERF_STOP(region)      /* This space intentionally left blank */ ;
#  define LCAPERF_PRINT             /* This space intentionally left blank */ ;

#endif

#endif /* PERFORMANCE_H */

