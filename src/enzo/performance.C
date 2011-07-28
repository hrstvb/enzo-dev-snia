/*****************************************************************************
 *                                                                           *
 * Copyright 2011 James Bordner                                              *
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


// list of regions and corresponding metric group

enum op_type {
  op_avg, // Average counter values over all processors
  op_sum, // Average counter values over all processors
  op_max, // Maximum counter value over all processors
  op_eff  // Efficiency (avg / max) of counter values over all procs
};
  
#define NUM_REGIONS 50
const struct {
  int          group;   // index of metric group;
  const char * name;    // name of region
} region_list[NUM_REGIONS] = 
  {
    {0, "EvolveLevel"},
    {0, "EvolveHierarchy"},
    {1, "ComovingExpansionTerms"},
    {1, "ComputePotentialFieldLevelZero"},
    {1, "ComputeRandomForcingNormalization"},
    {1, "CopyOverlappingMassField"},
    {1, "CreateFluxes"},
    {1, "DepositParticleMassField"},
    {1, "EvolvePhotons"},
    {1, "GetProjectedBoundaryFluxes"},
    {1, "grid_AddRandomForcing"},
    {1, "grid_ComputeAccelerationFieldExternal"},
    {1, "grid_CopyBaryonFieldToOldBaryonField"},
    {1, "grid_MultiSpeciesHandler"},
    {1, "grid_ParticleSplitter"},
    {1, "grid_SolveForPotential"},
    {1, "grid_SolveHydroEquations"},
    {1, "grid_StarParticleHandler"},
    {1, "InlineHaloFinder"},
    {1, "PrepareDensityField"},
    {1, "PrepareGravitatingMassField1"},
    {1, "PrepareGravitatingMassField2a"},
    {1, "PrepareGravitatingMassField2b"},
    {1, "ProjectSolutionToParentGrid"},
    {1, "RadiationFieldUpdate"},
    {1, "RebuildHierarchy"},
    {1, "SetBC_Parent"},
    {1, "SetBC_Siblings"},
    {1, "SetBoundaryConditions"},
    {1, "SetLevelTimeStep"},
    {1, "SolveForPotential"},
    {1, "star_FindFeedbackSphere"},
    {1, "star_FindFeedbackSphere2"},
    {1, "star_FindFeedbackSphere_Sum"},
    {1, "star_FindFeedbackSphere_Zero"},
    {1, "StarParticleAccretion"},
    {1, "StarParticleAddFeedback"},
    {1, "StarParticleDeath"},
    {1, "StarParticleFinalize"},
    {1, "StarParticleInitialize"},
    {1, "StarParticleSubtractAccretedMass"},
    {1, "star_SphereContained"},
    {1, "star_UpdatePositionVelocity"},
    {1, "UpdateFromFinerGrids"},
    {1, "UpdateParticlePositions"},
    {1, "UpdateStarParticleCount"}
  };

#define NUM_GROUPS   2 // maximum number of metric groups
#define NUM_METRICS 20 // maximum number of metrics per group

const struct {
  const char * name;    // name of derived metric for output file 
  const char * group;   // name of counter group of source counter
  const char * counter; // name of source counter
  const double scaling; // scaling factor
  const op_type op;
  
}  metric_list[NUM_GROUPS][NUM_METRICS] = { 
  // group 0  (AMR, time, flops, mpi, memory)
  {
    {"amr-zones-avg",       "user",  "count-zones",     1.0,  op_avg},
    {"amr-grids-avg",       "user",  "count-grids",     1.0,  op_avg},
    {"amr-ghosts-avg",      "user",  "count-ghosts",    1.0,  op_avg},
    {"amr-particles-avg",   "user",  "count-particles", 1.0,  op_avg},
    {"time-avg",            "basic", "time",            1e-6, op_avg},
    {"gflops-avg",          "papi",  "papi-fp-ops",     1e-9, op_avg},
    {"gflops-eff",          "papi",  "papi-fp-ops",     1.0,  op_eff},
    {"mpi-time-avg",        "mpi",   "mpi-time",        1e-6, op_avg},
    {"mpi-time-sync-avg",   "mpi",   "mpi-sync-time",   1e-6, op_avg},
    {"mpi-send-mbytes-avg", "mpi",   "mpi-send-bytes",  1e-6, op_avg},
    {"mpi-recv-mbytes-avg", "mpi",   "mpi-recv-bytes",  1e-6, op_avg},
    {"mpi-send-mbytes-eff", "mpi",   "mpi-send-bytes",  1.0,  op_eff},
    {"mpi-recv-mbytes-eff", "mpi",   "mpi-recv-bytes",  1.0,  op_eff},
    {"mem-curr-mbytes-avg", "mem",   "mem-curr-bytes",  1e-6, op_avg},
    {"mem-high-mbytes-avg", "mem",   "mem-high-bytes",  1e-6, op_avg},
    {"mem-curr-mbytes-eff", "mem",   "mem-curr-bytes",  1.0,  op_eff},
    {"mem-high-mbytes-eff", "mem",   "mem-high-bytes",  1.0,  op_eff}
  },
  // group 1 (time, flops, mpi)
  {
    {"time-avg",            "basic", "time",            1e-6, op_avg},
    {"gflops-avg",          "papi",  "papi-fp-ops",     1e-9, op_avg},
    {"gflops-eff",          "papi",  "papi-fp-ops",     1.0,  op_eff},
    {"mpi-send-mbytes-avg", "mpi",   "mpi-send-bytes",  1e-6, op_avg},
    {"mpi-recv-mbytes-avg", "mpi",   "mpi-recv-bytes",  1e-6, op_avg},
    {"mpi-send-mbytes-eff", "mpi",   "mpi-send-bytes",  1.0,  op_eff},
    {"mpi-recv-mbytes-eff", "mpi",   "mpi-recv-bytes",  1.0,  op_eff},
    {"mpi-time-avg",        "mpi",   "mpi-time",        1e-6, op_avg},
    {"mpi-time-sync-avg",   "mpi",   "mpi-sync-time",   1e-6, op_avg}
  }
};

FILE * fp_metric[NUM_REGIONS][NUM_METRICS] = {{0}};

void lcaperfInitialize (int max_level)
{

  // initialize lcaperf

  lcaperf.initialize ("out.lcaperf");

  // define lcaperf attributes

  lcaperf.new_attribute ("cycle", LCAP_INT);
  lcaperf.new_attribute ("level", LCAP_INT);

  // define lcaperf regions

  for (size_t i_region = 0; i_region<NUM_REGIONS; i_region++) {
    const char * region = region_list[i_region].name;
    if (region) lcaperf.new_region(region);
  }

  // define lcaperf counters

  lcaperf.new_counter ("count-zones",      counter_type_absolute);
  lcaperf.new_counter ("count-ghosts",     counter_type_absolute);
  lcaperf.new_counter ("count-grids",      counter_type_absolute);
  lcaperf.new_counter ("count-particles",  counter_type_absolute);

  // for (size_t level=0; level <= max_level; level++) {

  //   char lcaperf_counter_name[30];

  //   sprintf (lcaperf_counter_name,"count-zones-%"ISYM,level);
  //   lcaperf.new_counter(lcaperf_counter_name,counter_type_absolute);

  //   sprintf (lcaperf_counter_name,"count-ghosts-%"ISYM,level);
  //   lcaperf.new_counter(lcaperf_counter_name,counter_type_absolute);

  //   sprintf (lcaperf_counter_name,"count-grids-%"ISYM,level);
  //   lcaperf.new_counter(lcaperf_counter_name,counter_type_absolute);

  //   sprintf (lcaperf_counter_name,"count-particles-%"ISYM,level);
  //   lcaperf.new_counter(lcaperf_counter_name,counter_type_absolute);

  // }

  // create lcaperf directory and open fp_metric[][] files

  Eint32 ip = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&ip);
#endif

  if (ip == 0) {

    mkdir ("lcaperf",0777);
    chdir ("lcaperf");

    // open all performance metric files
    for (size_t i_region = 0; i_region < NUM_REGIONS; ++i_region) {

      int i_group         = region_list[i_region].group;
      const char * region = region_list[i_region].name;

      if (region) {
	for (size_t i_metric = 0; i_metric < NUM_METRICS; ++i_metric) {
	  const char * metric = metric_list[i_group][i_metric].name;
	  if (metric) {
	    char filename[80];
	    sprintf (filename,"%s.%s",region,metric);
	    fp_metric[i_region][i_metric] = fopen (filename,"w");
	  }
	}
      }
    }
    chdir ("..");
  }

  lcaperf.begin();
}

//----------------------------------------------------------------------

void lcaperfFinalize ()
{
  Eint32 ip = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&ip);
#endif

  if (ip == 0) {
    //  count regions for next loop
    for (size_t i_region = 0; i_region < NUM_REGIONS; ++i_region) {
      int i_group = region_list[i_region].group;
      const char * region = region_list[i_region].name;
      if (region) {
	for (size_t i_metric = 0; i_metric < NUM_METRICS; ++i_metric) {
	  const char * metric = metric_list[i_group][i_metric].name;
	  if (metric) {
	    fclose (fp_metric[i_region][i_metric]);
	  }
	}
      }
    }
  }
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

  size_t         cycle_index  = attributes_.index("cycle");
  std::string cycle_string = attributes_.value(cycle_index);

  // Allocate arrays for summing counters across levels if needed
  // and reducing along processors (based on metric_list[][].op)

  double value_reduce;
  long long counter_reduce_sum, counter_reduce_max;

  // LOOP OVER REGIONS

  for (size_t i_region = 0; i_region < NUM_REGIONS; ++i_region) {

    bool empty = true;
    
    const char * region = region_list[i_region].name;
    int         i_group = region_list[i_region].group;

    if (!region) continue;

    //    NOTE: EvolveLevel must be handled differently since it is recursive

    bool is_recursive = (strcmp(region,"EvolveLevel")==0);

    //    Only use level 0 for recursive functions since times are inclusive
    std::string level_string = (is_recursive) ? "0" : "*";

    // Create the region key to check counter keys against
    // [WARNING: DEPENDENCY ON NUMBER OF ATTRIBUTES AND ATTRIBUTE ORDERING]

    std::string region_key = 
      std::string(region) + ":" + cycle_string + ":" + level_string;

    // Loop over metrics

    for (size_t i_metric = 0; i_metric < NUM_METRICS; ++i_metric ) {

      const char * group   = metric_list[i_group][i_metric].group;
      const char * counter = metric_list[i_group][i_metric].counter;
      const double scaling = metric_list[i_group][i_metric].scaling;
      const op_type  op    = metric_list[i_group][i_metric].op;

      if (! group) continue;

      if (group && counters_.find(group) != counters_.end()) {

	// Loop over keys in the Counters object to sum counters over levels

	counter_reduce_sum = 0;
	counter_reduce_max = 0;

	ItCounterKeys itKeys (counters_[group]);
	size_t i_time = counters_[group]->index(counter);

	while (const char * key = ++itKeys) {

	  // Select matching keys

	  bool keys_match = attributes_.keys_match(key,region_key);

	  // Sum over levels
	  if (keys_match) {
	    long long * counters = itKeys.value();
	    counter_reduce_sum += counters[i_time];
	    counter_reduce_max += counters[i_time];
	  }
	}

	// Compute average and maximum over all processors

#ifdef USE_MPI
	switch (op) {
	case op_avg:
	case op_sum:
	  MPI_Allreduce (MPI_IN_PLACE,&counter_reduce_sum,1,
			 MPI_LONG_LONG,  MPI_SUM,  MPI_COMM_WORLD);
	  break;
	case op_max:
	  MPI_Allreduce (MPI_IN_PLACE,&counter_reduce_max,1,
			 MPI_LONG_LONG,  MPI_MAX,  MPI_COMM_WORLD);
	case op_eff:
	  MPI_Allreduce (MPI_IN_PLACE,&counter_reduce_sum,1,
			 MPI_LONG_LONG,  MPI_SUM,  MPI_COMM_WORLD);
	  MPI_Allreduce (MPI_IN_PLACE,&counter_reduce_max,1,
			 MPI_LONG_LONG,  MPI_MAX,  MPI_COMM_WORLD);
	  break;
	default:
	  // NOP
	  break;
	}
#endif

	// Set default reduction value

	value_reduce = (op==op_eff) ? (1.0) : (0.0);

	if (counter_reduce_sum != 0 || counter_reduce_max != 0) {

	  empty = false;

	  switch (op) {
	  case op_avg:
	    value_reduce = counter_reduce_sum/np;
	    break;
	  case op_sum:
	    value_reduce = counter_reduce_sum;
	    break;
	  case op_max:
	    value_reduce = counter_reduce_max;
	    break;
	  case op_eff:
	    value_reduce = (1.0*counter_reduce_sum/ np) / counter_reduce_max ;
	    break;
	  }
	  // Scale the value
	  value_reduce *= scaling;
	}

	// Print to file
	FILE * fp = fp_metric[i_region][i_metric];
	if (fp) {
	  fprintf (fp,"%lf",value_reduce);
	  fprintf (fp,"\n");
	  fflush(fp);
	}
      }
    }
  }

  // // Clear counters when done
  // for (size_t i_group=0; i_group<NUM_GROUPS; i_group++) {
  //   for (size_t i_metric = 0; i_metric<NUM_METRICS; ++i_metric) {
  //     const char * group = metric_list[i_group][i_metric].group;
  //     if (group && counters_[group]) counters_[group]->clear();
  //   }
  // }

}

#endif /* USE_LCAPERF */


