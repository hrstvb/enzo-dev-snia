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
// WARNING: macros_and_parameters.h may redefine "int" to "long int"
// WARNING: and "float" to "double" !
#include "macros_and_parameters.h"

#ifdef USE_LCAPERF

void lcaperfInitialize (int max_level)
{

  // Initialize lcaperf

  lcaperf.initialize ("out.lcaperf");

  // Define lcaperf attributes

  lcaperf.new_attribute ("cycle", LCAP_INT);
  lcaperf.new_attribute ("level", LCAP_INT);

  // Define lcaperf counters

  lcaperf.new_counter ("count-zones",      user_counter_type_absolute);
  lcaperf.new_counter ("count-ghosts",     user_counter_type_absolute);
  lcaperf.new_counter ("count-grids",      user_counter_type_absolute);
  lcaperf.new_counter ("count-particles",  user_counter_type_absolute);
  lcaperf.new_counter ("time-sim",         user_counter_type_absolute);

  // Select which regions to print()

  lcaperf.new_region ("EvolveHierarchy");
  lcaperf.new_region ("EvolveLevel");
  lcaperf.new_region ("ComovingExpansionTerms");
  lcaperf.new_region ("ComputePotentialFieldLevelZero");
  lcaperf.new_region ("ComputeRandomForcingNormalization");
  lcaperf.new_region ("CopyOverlappingMassField");
  lcaperf.new_region ("CreateFluxes");
  lcaperf.new_region ("DepositParticleMassField");
  lcaperf.new_region ("EvolvePhotons");
  lcaperf.new_region ("GetProjectedBoundaryFluxes");
  lcaperf.new_region ("grid_AddRandomForcing");
  lcaperf.new_region ("grid_ComputeAccelerationFieldExternal");
  lcaperf.new_region ("grid_CopyBaryonFieldToOldBaryonField");
  lcaperf.new_region ("grid_MultiSpeciesHandler");
  lcaperf.new_region ("grid_ParticleSplitter");
  lcaperf.new_region ("grid_SolveForPotential");
  lcaperf.new_region ("grid_SolveHydroEquations");
  lcaperf.new_region ("grid_StarParticleHandler");
  lcaperf.new_region ("InlineHaloFinder");
  lcaperf.new_region ("PrepareDensityField");
  lcaperf.new_region ("PrepareGravitatingMassField1");
  lcaperf.new_region ("PrepareGravitatingMassField2a");
  lcaperf.new_region ("PrepareGravitatingMassField2b");
  lcaperf.new_region ("ProjectSolutionToParentGrid");
  lcaperf.new_region ("RadiationFieldUpdate");
  lcaperf.new_region ("RebuildHierarchy");
  lcaperf.new_region ("SetBC_Parent");
  lcaperf.new_region ("SetBC_Siblings");
  lcaperf.new_region ("SetBoundaryConditions");
  lcaperf.new_region ("SetLevelTimeStep");
  lcaperf.new_region ("SolveForPotential");
  lcaperf.new_region ("star_FindFeedbackSphere");
  lcaperf.new_region ("star_FindFeedbackSphere2");
  lcaperf.new_region ("star_FindFeedbackSphere_Sum");
  lcaperf.new_region ("star_FindFeedbackSphere_Zero");
  lcaperf.new_region ("StarParticleAccretion");
  lcaperf.new_region ("StarParticleAddFeedback");
  lcaperf.new_region ("StarParticleDeath");
  lcaperf.new_region ("StarParticleFinalize");
  lcaperf.new_region ("StarParticleInitialize");
  lcaperf.new_region ("StarParticleSubtractAccretedMass");
  lcaperf.new_region ("star_SphereContained");
  lcaperf.new_region ("star_UpdatePositionVelocity");
  lcaperf.new_region ("UpdateFromFinerGrids");
  lcaperf.new_region ("UpdateParticlePositions");
  lcaperf.new_region ("UpdateStarParticleCount");

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
  lcaperf.end();
  lcaperf.finalize();
}

//======================================================================

#ifdef USE_MPI
const bool l_mpi = true;
#else
const bool l_mpi = false;
#endif

#ifdef USE_PAPI
const bool l_papi = true;
#else
const bool l_papi = false;
#endif

#ifdef USE_MEM
const bool l_mem = true;
#else
const bool l_mem = false;
#endif

#ifdef USE_HDF5
const bool l_disk = true;
#else
const bool l_disk = false;
#endif

//======================================================================
LcaPerfEnzo lcaperf;
//======================================================================

//----------------------------------------------------------------------
void LcaPerfEnzo::header ()
//----------------------------------------------------------------------
{
  if (fp_) {
    fprintf (fp_,"lcaperf: ");
    fprintf             (fp_,"   time(s)   " "   ");
    if (l_mpi)  fprintf (fp_,"   MPI(s)    " "   ");
    if (l_papi) fprintf (fp_," flops(GF)   " "   ");
    if (l_mem)  fprintf (fp_," memory(GB)  " "   ");
    if (l_disk) fprintf (fp_,"  disk(GB)   " "   ");
    fprintf (fp_,"cycle ");
    fprintf (fp_,"region\n");
  }
}

//----------------------------------------------------------------------
void LcaPerfEnzo::print ()
//----------------------------------------------------------------------
{
  TRACE("print");

  // Get comm size, required for computing average times etc. over all
  // processes
  Eint32 np;
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  int         cycle_index  = attributes_.index("cycle");
  std::string cycle_string = attributes_.value(cycle_index);

  // Allocate arrays for summing counters across levels and reducing
  // along processors (avg and max)

  const int i_avg = 0, i_max = 1;
  double time_avg, time_max, time_eff;

  long long counter_array_reduce[2];

  // Line to print
  std::string line;
  char field[80];

  // Loop over regions to print

  for (size_t i_region = 0; i_region < regions_.size(); ++i_region) {

    bool empty = true;
    
    if (fp_) line = "lcaperf: ";

    std::string region = regions_[i_region];

    // Create the region key to check counter keys against

    //    NOTE: EvolveLevel must be handled differently since it is recursive

    bool is_recursive = (region == "EvolveLevel");

    //    Only use level 0 for recursive functions since times are inclusive
    std::string level_string = (is_recursive) ? "0" : "*";

    //    WARNING: dependency on number of attributes and attribute ordering

    std::string region_key = region + ":" + cycle_string + ":" + level_string;

    //--------------------------------------------------
    // TOTAL TIME
    //--------------------------------------------------

    // Loop over keys in the Counters object to sum counters over levels

    counter_array_reduce[i_avg] = 0;
    counter_array_reduce[i_max] = 0;

    ItCounterKeys itKeys (counters_["basic"]);
    int i_time = counters_["basic"]->index("time");

    while (const char * key = ++itKeys) {

      // Select matching keys

      bool keys_match = attributes_.keys_match(key,region_key);

      if (keys_match) {
	// Get array of counters
	long long * counter_array = itKeys.value();
	// Sum over levels
	counter_array_reduce[i_avg] += counter_array[i_time];
	counter_array_reduce[i_max] += counter_array[i_time];
      }
    }

#ifdef USE_MPI
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_avg],1,MPI_LONG_LONG,
		   MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_max],1,MPI_LONG_LONG,
		   MPI_MAX,MPI_COMM_WORLD);
#endif

    time_avg = 0.0;
    time_eff = 1.0;

    if (counter_array_reduce[i_max] != 0) {

      empty = false;

      time_avg = 1e-6*counter_array_reduce[i_avg]/np;
      time_max = 1e-6*counter_array_reduce[i_max];
      time_eff = time_avg / time_max;
    }

    sprintf (field, "%06.2f %5.4f   ",time_avg,time_eff);

    line = line + field;

    //--------------------------------------------------
    // MPI TIME
    //--------------------------------------------------

    counter_array_reduce[i_avg] = 0;
    counter_array_reduce[i_max] = 0;

    ItCounterKeys itKeysMpi (counters_["mpi"]);
    i_time = counters_["mpi"]->index("mpi-time");

    while (const char * key = ++itKeysMpi) {

      // Select matching keys

      bool keys_match = attributes_.keys_match(key,region_key);

      if (keys_match) {
	// Get array of counters
	long long * counter_array = itKeysMpi.value();
	// Sum over levels
	counter_array_reduce[i_avg] += counter_array[i_time];
	counter_array_reduce[i_max] += counter_array[i_time];
      }
    }

#ifdef USE_MPI
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_avg],1,MPI_LONG_LONG,
		   MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_max],1,MPI_LONG_LONG,
		   MPI_MAX,MPI_COMM_WORLD);
#endif

    time_avg = 0.0;
    time_eff = 1.0;

    if (counter_array_reduce[i_max] != 0) {

      empty = false;

      time_avg = 1e-6*counter_array_reduce[i_avg]/np;
      time_max = 1e-6*counter_array_reduce[i_max];
      time_eff = time_max ? (time_avg / time_max) : 1.0;

    }

    sprintf (field, "%06.2f %5.4f   ",time_avg,time_eff);

    line = line + field;

    if (fp_ && ! empty) {

      sprintf (field , "%s   %s\n",cycle_string.c_str(),region.c_str());

      line = line + field;

      fprintf (fp_, line.c_str());
    }
  }

  //  counters_["basic"]->print(fp_ip_);
  //  counters_["mpi"]->print(fp_ip_);

  if (fp_) fprintf (fp_,"\n");

  counters_["basic"]->clear();
  counters_["mpi"]->clear();

  fflush(fp_);
}

#endif /* USE_LCAPERF */


