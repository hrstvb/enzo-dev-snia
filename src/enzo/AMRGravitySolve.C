#ifdef AMR_SOLVE

/// @file     AMRGravitySolve.C
/// @date     2010-12-29
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel Reynolds (reynolds@smu.edu)
/// @brief    Interface between Enzo and the AMRsolve linear solver
///
/// AMRGravitySolve serves as the single function interface between
/// Enzo and the AMRsolve linear solver package.  It solves the
/// Poisson equation on an Enzo AMR hierarchy between two levels specified
/// by the "level_coarse" and "level_fine" parameters.  An extra scaling
/// factor "f_scale" is passed, which can be 1.0.  
/// The right hand side and solution arrays are held in the Grid class, and
/// are accessed directly by the AMRsolve_Hierarchy::enzo_attach_grav() function.

#include "AMRsolve.h"

#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "ExternalBoundary.h"
#include "ProtoSubgrid.h"
#include "GridList.h"
#include "Grid.h"
#include "LevelHierarchy.h"
#include "Hierarchy.h"
//#include "debug.h"

const bool trace_amrsolve = true;
const bool debug_amrsolve = true;

void AMRGravitySolve 
(
 LevelHierarchyEntry * LevelArray[],    // Enzo's current LevelArray
 double                f_scale,         // scaling for the right-hand side
 int                   level_coarse,    // coarsest level
 int                   level_fine       // finest level
 )
{

  if (trace_amrsolve && pmpi->is_root()) {
    printf("%s:%d %d AMRGravitySolve(%d %d)\n",
	   __FILE__,__LINE__,pmpi->ip(),
	   Eint32(level_coarse),Eint32(level_fine));
    fflush(stdout);
  }

  LCAPERF_START("amr_solve");
 
  // Insert Enzo grids in this level into a AMRsolve hierarchy
  AMRsolve_Hierarchy* hierarchy = new AMRsolve_Hierarchy;

  LCAPERF_START("amrsolve_attach_grav");
  hierarchy->enzo_attach_grav(LevelArray,level_coarse,level_fine);
  LCAPERF_STOP("amrsolve_attach_grav");

  // Initialize the AMRsolve hierarchy
  AMRsolve_Domain domain(3, DomainLeftEdge, DomainRightEdge);
  bool is_periodic = true;

  LCAPERF_START("amrsolve_hierarchy");
  hierarchy->initialize(domain,*pmpi,is_periodic);
  LCAPERF_STOP("amrsolvehierarchy");

  // Initialize the AMRsolve linear system
  LCAPERF_START("amrsolve_matrix");
  AMRsolve_Hypre_Grav amrsolve(*hierarchy, *amrsolve_parameters);
  amrsolve.init_hierarchy(*pmpi);
  amrsolve.init_stencil();
  amrsolve.init_graph();
  std::vector<AMRsolve_Point *> points; // ignored
  amrsolve.init_elements(points, f_scale);
  LCAPERF_STOP("amrsolve_matrix");

  // Solve the linear system
  LCAPERF_START("amrsolve_solve");
  amrsolve.solve();
  LCAPERF_STOP("amrsolve_solve");

  // Display solver results
  amrsolve.evaluate();

  // Output
  static int cycle = 0;

  if (cycle % 10 == 0) {
    char filename [80];

    ItHierarchyGridsLocal itgl (*hierarchy);
    printf("%s:%d scale = %g\n",__FILE__,__LINE__,f_scale);

    while (AMRsolve_Grid* grid = itgl++) {
      sprintf(filename,"XH.id=%d.proc=%ld.cycle=%d",
	      grid->id(),MyProcessorNumber,cycle);
      printf("%s:%d writing %s\n",__FILE__,__LINE__,filename);
      FILE* fp = fopen(filename,"w");
      Eint32 nx,ny,nz;
      Scalar* values = grid->get_u(&nx,&ny,&nz);
      for (int ix=0; ix<nx; ix++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int iz=0; iz<nz; iz++) {
	    int i = ix + nx * ( iy + ny * iz);
	    fprintf(fp,"%d %d %d %g\n",ix,iy,iz,values[i]);
	  }
	}
      }
      fclose(fp);
    }
  }
  cycle ++;

  // Clean up
  hierarchy->enzo_detach();

  LCAPERF_STOP("amr_solve");

  delete hierarchy;
  hierarchy = NULL;
}

#endif   // AMR_SOLVE
