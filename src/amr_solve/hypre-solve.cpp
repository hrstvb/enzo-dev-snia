/// @file      hypre-solve.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Main driver for hypre AMR test solver

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include <map>
#include <string>
#include <vector>


#include "HYPRE_sstruct_ls.h"

#include "AMRsolve_defs.h"

//----------------------------------------------------------------------

#define debug true
#define trace false

//----------------------------------------------------------------------

#include "AMRsolve_scalar.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_error.h"
#include "AMRsolve_performance.h"
#include "AMRsolve_point.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_problem.h"
#include "AMRsolve_hypre.h"

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  // Determine executable name

  std::string exec_name (argv[0]);
  int size = exec_name.rfind("/");
  exec_name.replace(0,size+1,"");

  // --------------------------------------------------
  // MPI initialization
  // --------------------------------------------------

  if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);

  pmpi = new AMRsolve_Mpi(&argc,&argv);

  if (debug) printf ("DEBUG %s:%d mpi (ip,np) = (%d,%d)\n",
		     __FILE__,__LINE__,pmpi->ip(),pmpi->np());

  AMRsolve_Grid::set_mpi (*pmpi);


  if (argc==2) {

    // --------------------------------------------------
    // LCAPERF initialization
    // --------------------------------------------------

    LCAPERF_BEGIN("EL");

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    LCAPERF_START("problem");
    AMRsolve_Problem problem;
    LCAPERF_STOP("problem");

    LCAPERF_START("problem-read");
    problem.read  (filename);
    LCAPERF_STOP("problem-read");

    // --------------------------------------------------
    // Initialize the hierarchy
    // --------------------------------------------------

    AMRsolve_Hierarchy& hierarchy = problem.hierarchy();
    AMRsolve_Grid::set_domain(problem.domain());  // only needed by geomview viz
    AMRsolve_Level::set_domain(problem.domain()); // only needed by geomview viz

    // determine interconnections between grids

    bool is_periodic = problem.parameters().value("boundary") == "periodic";
    LCAPERF_START("hierarchy-initialize");
    hierarchy.initialize(problem.domain(), *pmpi,is_periodic);
    LCAPERF_STOP("hierarchy-initialize");

    if (debug) problem.print ();

    // WARNING: problem-size returns the size of the root-grid, not the entire hierarchy

    // --------------------------------------------------
    // Initialize hypre
    // --------------------------------------------------

    AMRsolve_Hypre hypre(hierarchy,problem.parameters());


    LCAPERF_START("hypre-init-hierarchy");
    hypre.init_hierarchy (*pmpi);
    LCAPERF_STOP("hypre-init-hierarchy");

    // Initialize the stencils
    
    LCAPERF_START("hypre-init-stencil");
    hypre.init_stencil ();
    LCAPERF_STOP("hypre-init-stencil");

    // Initialize the graph

    LCAPERF_START("hypre-init-graph");
    hypre.init_graph ();
    LCAPERF_STOP("hypre-init-graph");

    // Initialize the elements of matrix A and vector B

    LCAPERF_START("hypre-init-elements");
    hypre.init_elements (problem.points());
    LCAPERF_STOP("hypre-init-elements");

    // --------------------------------------------------
    // Solve the linear system A X = B
    // --------------------------------------------------

    LCAPERF_START("hypre-solve");
    hypre.solve ();
    LCAPERF_STOP("hypre-solve");

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    LCAPERF_START("hypre-evaluate");
    hypre.evaluate ();
    LCAPERF_STOP("hypre-evaluate");

    // --------------------------------------------------
    // lcaperf Finalize
    // --------------------------------------------------

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------


    MPI_Finalize();
    delete pmpi;

    LCAPERF_END("EL");

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

}


