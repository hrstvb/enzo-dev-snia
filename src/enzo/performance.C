/*****************************************************************************
 *                                                                           *
 * Copyright 2011 James Bordner
 * Copyright 2011 Laboratory for Computational Astrophysics                  *
 * Copyright 2011 Board of Trustees of the University of Illinois            *
 * Copyright 2011 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/

#include "performance.h"
// WARNING: macros_and_parameters.h may redefine "int" and "float" !
#include "macros_and_parameters.h"

#ifdef USE_LCAPERF

// LIST OF LCAPERF REGIONS TO OUTPUT
const char * region_list[] = {
  "EvolveHierarchy",
  "EvolveLevel",
  "ComovingExpansionTerms",
  "ComputePotentialFieldLevelZero",
  "ComputeRandomForcingNormalization",
  "CopyOverlappingMassField",
  "CreateFluxes",
  "DepositParticleMassField",
  "EvolvePhotons",
  "GetProjectedBoundaryFluxes",
  "grid_AddRandomForcing",
  "grid_ComputeAccelerationFieldExternal",
  "grid_CopyBaryonFieldToOldBaryonField",
  "grid_MultiSpeciesHandler",
  "grid_ParticleSplitter",
  "grid_SolveForPotential",
  "grid_SolveHydroEquations",
  "grid_StarParticleHandler",
  "InlineHaloFinder",
  "PrepareDensityField",
  "PrepareGravitatingMassField1",
  "PrepareGravitatingMassField2a",
  "PrepareGravitatingMassField2b",
  "ProjectSolutionToParentGrid",
  "RadiationFieldUpdate",
  "RebuildHierarchy",
  "SetBC_Parent",
  "SetBC_Siblings",
  "SetBoundaryConditions",
  "SetLevelTimeStep",
  "SolveForPotential",
  "star_FindFeedbackSphere",
  "star_FindFeedbackSphere2",
  "star_FindFeedbackSphere_Sum",
  "star_FindFeedbackSphere_Zero",
  "StarParticleAccretion",
  "StarParticleAddFeedback",
  "StarParticleDeath",
  "StarParticleFinalize",
  "StarParticleInitialize",
  "StarParticleSubtractAccretedMass",
  "star_SphereContained",
  "star_UpdatePositionVelocity",
  "UpdateFromFinerGrids",
  "UpdateParticlePositions",
  "UpdateStarParticleCount",
  // NULL STRING SIGNALS ARRAY END: REQUIRED
  ""
};

// LIST OF LCAPERF METRICS TO OUTPUT

const char * metric_list[] = { 
  // "counter group", "metric"
  "basic", "time",
  "mpi",   "mpi-time",
  "mpi",   "mpi-sync-time",
  "mpi",   "mpi-send-bytes",
  "mpi",   "mpi-recv-bytes",
  // NULL STRINGS SIGNALS ARRAY END: REQUIRED
  "","" 
};

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

  for (int i=0; strlen(region_list[i])>0; i++) {
    lcaperf.new_region (region_list[i]);
  }

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


//======================================================================
LcaPerfEnzo lcaperf;
//======================================================================

//----------------------------------------------------------------------
void LcaPerfEnzo::header ()
//----------------------------------------------------------------------
{
  const bool l_mpi  = counters_.find("mpi")  != counters_.end();
  const bool l_papi = counters_.find("papi") != counters_.end();
  const bool l_mem  = counters_.find("mem")  != counters_.end();
  const bool l_disk = counters_.find("disk") != counters_.end();

  if (fp_) {
    fprintf (fp_,"lcaperf: ");
    fprintf             (fp_,"    time(s)    " "   ");
    if (l_mpi) {
      fprintf (fp_,"   MPI time(s)    " "  ");
      fprintf (fp_," MPI time sync(s) " "  ");
      fprintf (fp_,"MPI send bytes(MB)" "  ");
      fprintf (fp_,"MPI recv bytes(MB)" "  ");
    }
    if (l_papi) fprintf (fp_,"  flops(GF)    " "   ");
    if (l_mem)  fprintf (fp_,"  memory(GB)   " "   ");
    if (l_disk) fprintf (fp_,"   disk(GB)    " "   ");
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
  Eint32 np = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&np);
#endif

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

  // LOOP OVER REGIONS

  for (size_t i_region = 0; i_region < regions_.size(); ++i_region) {

    bool empty = true;
    
    if (fp_) line = "lcaperf: ";

    std::string region = regions_[i_region];

    //    NOTE: EvolveLevel must be handled differently since it is recursive

    bool is_recursive = (region == "EvolveLevel");

    //    Only use level 0 for recursive functions since times are inclusive
    std::string level_string = (is_recursive) ? "0" : "*";

    // Create the region key to check counter keys against
    // [WARNING: dependency on number of attributes and attribute ordering]

    std::string region_key = region + ":" + cycle_string + ":" + level_string;

    
    // LOOP OVER METRICS

    for (int i_metric = 0; strlen(metric_list[i_metric]) > 0; /* EMPTY */ ) {

      const char * group  = metric_list[i_metric++];
      const char * metric = metric_list[i_metric++];

      // Loop over keys in the Counters object to sum counters over levels

      counter_array_reduce[i_avg] = 0;
      counter_array_reduce[i_max] = 0;

      ItCounterKeys itKeys (counters_[group]);
      int i_time = counters_[group]->index(metric);

      while (const char * key = ++itKeys) {

	// Select matching keys

	bool keys_match = attributes_.keys_match(key,region_key);

	// Sum over levels
	if (keys_match) {
	  long long * counter_array = itKeys.value();
	  counter_array_reduce[i_avg] += counter_array[i_time];
	  counter_array_reduce[i_max] += counter_array[i_time];
	}
      }

      // Compute average and maximum over all processors

#ifdef USE_MPI
      MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_avg],1,MPI_LONG_LONG,
		     MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_max],1,MPI_LONG_LONG,
		     MPI_MAX,MPI_COMM_WORLD);
#endif

      // 
      time_avg = 0.0;
      time_eff = 1.0;

      if (counter_array_reduce[i_max] != 0) {

	empty = false;

	time_avg = 1e-6*counter_array_reduce[i_avg]/np;
	time_max = 1e-6*counter_array_reduce[i_max];
	time_eff = time_avg / time_max;
      }

      sprintf (field, "%06.4f %5.6f     ",time_avg,time_eff);

      line = line + field;

    }

    // Print the line

    if (fp_ && ! empty) {
      fprintf (fp_ , "%s %s   %s\n",
	       line.c_str(),
	       cycle_string.c_str(),
	       region.c_str());
    }

  }

  if (fp_) fprintf (fp_,"\n");

  if (counters_.find("basic") != counters_.end()) counters_["basic"]->clear();
  if (counters_.find("mpi") != counters_.end())   counters_["mpi"]->clear();

  fflush(fp_);
}

#endif /* USE_LCAPERF */


