#ifdef HYPRE_GRAV

/// @file     HypreGravitySolve.C
/// @date     2009-09-02
/// @author   James Bordner (jobordner@ucsd.edu)
/// @brief    Interface between Enzo and the hypre-solve linear solver
///
/// HypreGravitySolve serves as the single function interface between
/// Enzo and the hypre-solve linear solver package.  It solves the
/// Poisson equation on an Enzo AMR hierarchy between two levels specified
/// by the "level_coarse" and "level_fine" parameters.  An extra scaling
/// factor "f_scale" is passed, which can be 1.0.  
/// The right hand side is specified by "Grid::hypre_grav_b", which is
/// created and initialized in the "Hierarchy::enzo_attach()" function
/// in [newgrav-]hierarchy.C

#include "newgrav.h"

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
#include "debug.h"

const bool trace_hypre = true;
const bool debug_hypre = true;

void HypreGravitySolve 
(
 LevelHierarchyEntry * LevelArray[],    // Enzo's current LevelArray
 double                f_scale,         // scaling for the right-hand side
 int                   level_coarse,    // coarsest level
 int                   level_fine       // finest level
 )
{

  if (trace_hypre && pmpi->is_root()) {
    printf ("%s:%d %d HypreGravitySolve(%d %d)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    Eint32(level_coarse),Eint32(level_fine));
    fflush(stdout);
  }

  LCAPERF_START("hypre_solve");
 
  // Insert Enzo grids in this level into a hypre-solve hierarchy

  Hierarchy * hierarchy = new Hierarchy;

  LCAPERF_START("hypre_attach");
  hierarchy->enzo_attach(LevelArray,level_coarse,level_fine);
  LCAPERF_STOP("hypre_attach");

  // Initialize the hypre-solve hierarchy

  Domain domain (3, DomainLeftEdge, DomainRightEdge);
  bool   is_periodic = true;

  LCAPERF_START("hypre_hierarchy");
  hierarchy->initialize(domain,*pmpi,is_periodic);
  LCAPERF_STOP("hypre_hierarchy");

  // Initialize the hypre-solve linear system

  LCAPERF_START("hypre_matrix");
  Hypre hypre (*hierarchy,*hypre_parameters);

  hypre.init_hierarchy (*pmpi);
  hypre.init_stencil ();
  hypre.init_graph ();
  std::vector<Point *>  points; // ignored
  hypre.init_elements (points,f_scale);
  LCAPERF_STOP("hypre_matrix");

  // Solve the linear system

  LCAPERF_START("hypre_solve");
  hypre.solve ();
  LCAPERF_STOP("hypre_solve");

  // Display solver results

  hypre.evaluate ();

  // Output

  static int cycle = 0;

  if (cycle % 10 == 0) {

    char filename [80];

    ItHierarchyGridsLocal itgl (*hierarchy);
    printf ("%s:%d scale = %g\n",__FILE__,__LINE__,f_scale);
    while (Grid * grid = itgl++) {
      sprintf (filename,"XH.id=%d.proc=%ld.cycle=%d",
	       grid->id(),MyProcessorNumber,cycle);
      printf ("%s:%d writing %s\n",__FILE__,__LINE__,filename);
      FILE * fp = fopen (filename,"w");
      Eint32 nx,ny,nz;
      Scalar * values = grid->get_u(&nx,&ny,&nz);
      for (int ix=0; ix<nx; ix++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int iz=0; iz<nz; iz++) {
	    int i = ix + nx * ( iy + ny * iz);
	    fprintf (fp,"%d %d %d %g\n",ix,iy,iz,values[i]);
	    
	  }
	}
      }
      fclose(fp);
    }
  }

  cycle ++;

  // Clean up

  hierarchy->enzo_detach();

  LCAPERF_STOP("hypre_solve");

  delete hierarchy;
  hierarchy = NULL;
}
#endif
