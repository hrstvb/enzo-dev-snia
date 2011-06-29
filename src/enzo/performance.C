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
//======================================================================
//
// File:        performance.C
//
// Description: Performance-related code
//
//----------------------------------------------------------------------
//
// Namespaces:  jb
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

#include "performance.h"
#include "macros_and_parameters.h"

#ifdef USE_LCAPERF

void lcaperfInitialize (int max_level)
{

  // Initialize lcaperf

  // Define lcaperf attributes

  lcaperf.new_attribute ("cycle", LCAP_INT);
  lcaperf.new_attribute ("level", LCAP_INT);

  // Define lcaperf counters

  lcaperf.new_counter ("count-zones",      user_counter_type_absolute);
  lcaperf.new_counter ("count-ghosts",     user_counter_type_absolute);
  lcaperf.new_counter ("count-grids",      user_counter_type_absolute);
  lcaperf.new_counter ("count-particles",  user_counter_type_absolute);
  lcaperf.new_counter ("time-sim",         user_counter_type_absolute);

  for (int level=0; level <= max_level; level++) {

    char lcaperf_counter_name[30];

    sprintf (lcaperf_counter_name,"count-zones-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-ghosts-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-grids-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

    sprintf (lcaperf_counter_name,"count-particles-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,user_counter_type_absolute);

  }

  lcaperf.begin();
}

//----------------------------------------------------------------------

void lcaperfFinalize ()
{
  lcaperf.print();
  lcaperf.end();
}

#endif /* USE_LCAPERF */


