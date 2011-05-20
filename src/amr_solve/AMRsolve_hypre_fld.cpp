/// @file      AMRsolve_hypre_fld.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Hypre_FLD class

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"
#include "HYPRE_krylov.h"

#include "hdf5.h"

#include "AMRsolve_defs.h"

//----------------------------------------------------------------------

// defines

#define WHERE printf ("%s:%d ",__FILE__,__LINE__);

//----------------------------------------------------------------------

// Constants

const bool debug = false;
const bool trace = false;

//----------------------------------------------------------------------

// Typedefs

typedef int int3[3];

//----------------------------------------------------------------------

#include "AMRsolve_mpi.h"
#include "AMRsolve_scalar.h"
#include "AMRsolve_error.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_hypre_fld.h"
#include "AMRsolve_error.h"

//======================================================================

// Coefficient for Poisson problem

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// AMRsolve_Hypre_FLD constructor
AMRsolve_Hypre_FLD::AMRsolve_Hypre_FLD(AMRsolve_Hierarchy& hierarchy, 
				       AMRsolve_Parameters& parameters)
  : grid_(0), graph_(0), stencil_(0), A_(0), B_(0), X_(0),
    solver_(0), parameters_(&parameters), hierarchy_(&hierarchy),
    resid_(-1.0), iter_(-1), r_factor_(const_r_factor), Nchem_(-1),
    theta_(-1.0), dt_(-1.0), aval_(-1.0), aval0_(-1.0), adot_(-1.0), 
    adot0_(-1.0), HIconst_(-1.0), HeIconst_(-1.0), HeIIconst_(-1.0), 
    nUn_(-1.0), nUn0_(-1.0), lUn_(-1.0), lUn0_(-1.0), rUn_(-1.0), rUn0_(-1.0)
{
  // set array-valued items
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = -1;
}

//----------------------------------------------------------------------

/// AMRsolve_Hypre_FLD destructor
AMRsolve_Hypre_FLD::~AMRsolve_Hypre_FLD()
{
  // destroy HYPRE objects that we created along the way
  char error_message[100];
  if (HYPRE_SStructVectorDestroy(B_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy B_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructVectorDestroy(X_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy X_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructMatrixDestroy(A_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy A_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructGraphDestroy(graph_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy graph_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructStencilDestroy(stencil_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy stencil_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructGridDestroy(grid_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy grid_\n");
    ERROR(error_message);
  }

}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy
/** Creates a hypre grid, with one part per level and one box per 
    Grid patch object, for an AMR problem.  Sets grid box extents, 
    grid part variables, and periodicity. */
void AMRsolve_Hypre_FLD::init_hierarchy(AMRsolve_Mpi& mpi)
{

  int dim       = hierarchy_->dimension();
  int num_parts = hierarchy_->num_levels();

  // Create the hypre grid
  _TRACE_;
  HYPRE_SStructGridCreate(MPI_COMM_WORLD, dim, num_parts, &grid_);

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level* level = itl++) {

    _TRACE_;
    int part = level->index();

    // Set extents for boxes that comprise the hypre grid
    ItLevelGridsLocal itgl (*level);
    while (AMRsolve_Grid* grid = itgl++) {
      int lower[3] = {grid->index_lower(0),
		      grid->index_lower(1),
		      grid->index_lower(2)};
      int upper[3] = {grid->index_upper(0),
		      grid->index_upper(1),
		      grid->index_upper(2)};
      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
    } // while grid = itgl++

    _TRACE_;
    // Create a single cell-centered variable for each grid part (level)
    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;
    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);

    // Set periodicity of the grid part
    int period[3] = { hierarchy_->period_index(0,part),
		      hierarchy_->period_index(1,part),
		      hierarchy_->period_index(2,part) };

    _TRACE_;
    HYPRE_SStructGridSetPeriodic(grid_, part, period);

  } // while level = itl++

  // When finished, assemble the hypre grid
  _TRACE_;

  HYPRE_SStructGridAssemble(grid_);
  _TRACE_;
  
} // AMRsolve_Hypre_FLD::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  
/** Creates and initializes a hypre stencil object. */
void AMRsolve_Hypre_FLD::init_stencil()
{

  _TRACE_;
  int dim = hierarchy_->dimension();

  _TRACE_;
  HYPRE_SStructStencilCreate(dim,dim*2+1,&stencil_);

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  _TRACE_;
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 0, entries[0], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 1, entries[1], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 2, entries[2], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry(stencil_, 3, entries[3], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry(stencil_, 4, entries[4], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry(stencil_, 5, entries[5], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry(stencil_, 6, entries[6], 0);
  _TRACE_;

} // AMRsolve_Hypre_FLD::init_stencil()

//----------------------------------------------------------------------

/// Initialize the graph.
/** Creates a graph containing the matrix nonzero structure.  Graph
    edges include both those for nonzeros from the stencil within each
    part (level), and nonzeros for graph entries connecting linked
    parts.  The matrix nonzero structure is generally nonsymmetric.
    Only the stencil step is required for unigrid problems. */
void AMRsolve_Hypre_FLD::init_graph()
{
  // Create the hypre graph object
  _TRACE_;
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid_, &graph_);
  
  _TRACE_;
  HYPRE_SStructGraphSetObjectType(graph_, HYPRE_SSTRUCT);

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level* level = itl++) {

    int part = level->index();

    // Define stencil connections within each level
    _TRACE_;
    HYPRE_SStructGraphSetStencil(graph_, part, 0, stencil_);

    // Define graph connections between levels
    _TRACE_;
    if (part > 0) {
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid* grid = itag++)  init_graph_nonstencil_(*grid);
    } // if part > 0

    // Initialize face counters for subsequent matrix inter-level entries
    _TRACE_;
    ItLevelGridsAll itag (*level);
    while (AMRsolve_Grid * grid = itag++) {
      int dim = hierarchy_->dimension();
      grid->init_counter(dim*2+1);
    } // while grid = itag++

  } // while level = itl++

  // Assemble the hypre graph
  _TRACE_;
  HYPRE_SStructGraphAssemble(graph_);
  _TRACE_;

} // AMRsolve_Hypre_FLD::init_graph()

//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b
/* Creates a matrix with a given nonzero structure, and sets nonzero
   values. */
void AMRsolve_Hypre_FLD::init_elements(double dt, int Nchem, double theta, 
				       double aval, double aval0, 
				       double adot, double adot0, 
				       double HIconst, double HeIconst, 
				       double HeIIconst, double nUn, 
				       double nUn0, double lUn, 
				       double lUn0, double rUn, 
				       double rUn0, int BdryType[3][2])
{
  // set input arguments into AMRsolve_Hypre_FLD object
  dt_        = dt;
  Nchem_     = Nchem;
  theta_     = theta;
  aval_      = aval;
  aval0_     = aval0;
  adot_      = adot;
  adot0_     = adot0;
  HIconst_   = HIconst;
  HeIconst_  = HeIconst;
  HeIIconst_ = HeIIconst;
  nUn_       = nUn;
  nUn0_      = nUn0;
  lUn_       = lUn;
  lUn0_      = lUn0;
  rUn_       = rUn;
  rUn0_      = rUn0;
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = BdryType[i][j];


  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph_, &A_);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &X_);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &B_);

  // Set the object types
  if (parameters_->value("solver") == "bicgstab-boomer") {
    HYPRE_SStructMatrixSetObjectType(A_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType(X_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType(B_,HYPRE_PARCSR);
  } else {
    HYPRE_SStructMatrixSetObjectType(A_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType(X_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType(B_,HYPRE_SSTRUCT);
  }

  // Initialize the hypre matrix and vector objects
  HYPRE_SStructMatrixInitialize(A_);
  HYPRE_SStructVectorInitialize(X_);
  HYPRE_SStructVectorInitialize(B_);

  //--------------------------------------------------
  // Initialize the matrix A_ elements
  //--------------------------------------------------
  init_elements_matrix_();

  //--------------------------------------------------
  // Initialize B_ elements 
  //--------------------------------------------------
  init_elements_rhs_();

  // Assemble the matrix and vectors
  HYPRE_SStructMatrixAssemble(A_);
  HYPRE_SStructVectorAssemble(B_);
  HYPRE_SStructVectorAssemble(X_);

  // Optionally write the matrix and vectors to a file for debugging
  if (parameters_->value("dump_a") == "true")  HYPRE_SStructMatrixPrint("A-hypre",A_,0);
  if (parameters_->value("dump_x") == "true")  HYPRE_SStructVectorPrint("X0-hypre",X_,0);
  if (parameters_->value("dump_b") == "true")  HYPRE_SStructVectorPrint("B-hypre",B_,0);

} // AMRsolve_Hypre_FLD::init_elements()

//----------------------------------------------------------------------

/// Initialize and solve the linear solver
void AMRsolve_Hypre_FLD::solve()
{
  std::string solver = parameters_->value("solver");
  int         levels = hierarchy_->num_levels();

  int    itmax  = 0;
  double restol = 0.0;

  // Check solver parameters
  std::string sitmax  = parameters_->value("solver_itmax");
  std::string srestol = parameters_->value("solver_restol");

  // If not defined, then define them
  if (sitmax == "")  parameters_->add_parameter("solver_itmax","200");
  if (srestol == "") parameters_->add_parameter("solver_restol","1e-6");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  if (solver == "pfmg" && levels == 1) {
    solve_pfmg_(itmax,restol);

  } else if (solver == "fac"  && levels > 1) {
    solve_fac_(itmax,restol);

  } else if (solver == "bicgstab") {
    solve_bicgstab_(itmax,restol);

  } else if (solver == "bicgstab-boomer") {
    solve_bicgstab_boomer_(itmax,restol);

  } else if (solver == "gmres") {
    solve_gmres_(itmax,restol);

  } else {
    char error_message[100];
    sprintf(error_message, "AMRsolve_Hypre_FLD::solve called with illegal "
	    "combination of solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters_->value("dump_x") == "true")  
    HYPRE_SStructVectorPrint("X-hypre",X_,1);

} // AMRsolve_Hypre_FLD::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve, return values (0=success, 1=failure)
int AMRsolve_Hypre_FLD::evaluate()
{
  // check whether solution/rhs was requested
  if (parameters_->value("dump_x") == "true" || 
      parameters_->value("dump_b") == "true") {

    // iterate over processor-local grids
    ItHierarchyGridsLocal itg(*hierarchy_);
    while (AMRsolve_Grid* grid = itg++) {

      char filename[80];

      // get level & grid information
      int level = grid->level();
      int lower[3],upper[3];
      grid->get_limits(lower,upper);
      
      if (parameters_->value("dump_x") == "true") {
	// extract Enzo solution
	int nx[3];
	HYPRE_SStructVectorGather(X_);
	HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0,
					grid->get_u(&nx[0],&nx[1],&nx[2]));  
	sprintf(filename,"X.%d",grid->id());
	grid->write("header",filename);
	grid->write("u",filename);
      }
      
      if (parameters_->value("dump_b") == "true") {
	// extract Enzo rhs
	int nb[3];
	HYPRE_SStructVectorGather(B_);
	HYPRE_SStructVectorGetBoxValues(B_, level, lower, upper, 0,
					grid->get_f(&nb[0],&nb[1],&nb[2]));  
	sprintf(filename,"B.%d",grid->id());
	grid->write("f",filename);
      }
      
    } // grid = itg++
  } // if dump_x or dump_b


  // check for error flags; output info to stdout, if requested
  int err_flag = 0;

  // Residual too high
  double restol = atof(parameters_->value("solver_restol").c_str());
  if (resid_ > restol) {
    if (pmpi->is_root()) 
      printf("Diverged: %g > %g\n", resid_,restol);
    err_flag = 1;
  }

  // Iterations reached limit
  int itmax = atoi(parameters_->value("solver_itmax").c_str());
  if (iter_ >= itmax) {
    if (pmpi->is_root()) 
      printf("Stalled: %d >= %d\n", iter_,itmax);
    err_flag = 1;
  }

//   // Appears to have completed successfully
//   if (err_flag == 0)  
//     if (pmpi->is_root()) {
//       printf("AMRsolve_Hypre_FLD Success!\n"); 
//       fflush(stdout); 
//     }

  return err_flag;

} // AMRsolve_Hypre_FLD::evaluate()

//----------------------------------------------------------------------

/// Extracts HYPRE solution and updates Enzo radiation field
void AMRsolve_Hypre_FLD::update_enzo()
{
  // iterate over grids on this processor
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // get level, grid information
    int level = grid->level();
    int lower[3],upper[3];
    grid->get_limits(lower,upper);

    // extract Enzo solution
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    HYPRE_SStructVectorGather(X_);
    HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0, u);  

    // access Enzo radiation field
    Scalar* E = grid->get_E();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_Ghosts(ghosts);
    int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];
    int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];
    int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];

    // update Enzo data with amrsolve solution
    int k0, k1, k2, k, i0, i1, i2, i;
    for (i2=0; i2<n3[2]; i2++) {
      k2 = ghosts[2][0] + i2;

      for (i1=0; i1<n3[1]; i1++) {
	k1 = ghosts[1][0] + i1;

	for (i0=0; i0<n3[0]; i0++) {
	  k0 = ghosts[0][0] + i0;

	  // compute indices of amrsolve, enzo cells
	  i = i0 + n3[0]*(i1 + n3[1]*i2);
	  k = k0 + en0*(k1 + en1*k2);

	  // update Enzo solution with amrsolve solution (corrector)
	  E[k] += u[i];
	}
      }
    }

  } // grid = itg++

} // AMRsolve_Hypre_FLD::update_enzo()


//----------------------------------------------------------------------

/// dumps HYPRE matrix and RHS (called when aborting solve)
void AMRsolve_Hypre_FLD::abort_dump()
{

  // have HYPRE dump out everything it knows to disk
  HYPRE_SStructMatrixPrint("A.mat",A_,0);
  HYPRE_SStructVectorPrint("x.vec",X_,0);
  HYPRE_SStructVectorPrint("b.vec",B_,0);

} // AMRsolve_Hypre_FLD::abort_dump()


//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_() with phase == phase_graph to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_() with phase
/// == phase_matrix to set nonstencil matrix entries.
void AMRsolve_Hypre_FLD::init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase)
{
  // Input parameter check
  if ( !(phase == phase_graph || phase == phase_matrix) ) {
    char error_message[80];
    sprintf(error_message,"init_matrix_nonstencil_ called with phase = %d",
	    int(phase));
    ERROR(error_message);
  } // if phase unexpected

  // global grid index limits
  int index_global[3][2];
  grid.indices(index_global);
  
  // Determine the discretization method (constant, linear, quadratic)
  // ***Currently only constant is implemented***
  enum discret_type_enum {discret_type_unknown, discret_type_const};
  discret_type_enum discret_type;

  if (parameters_->value("discret") == "constant") 
    discret_type = discret_type_const;

  // Loop over each face zone in the grid, adding non-stencil graph
  // entries wherever a zone is adjacent to a coarse zone.  Both
  // fine-to-coarse and coarse-to-fine entries are added.
  int level_fine   = grid.level();
  int level_coarse = grid.level() - 1;
  assert (level_coarse >= 0);

  for (int axis=0; axis<3; axis++) {

    // axis0:       axis normal to face
    // axis1,axis2: axes within face
    int axis0 = axis;
    int axis1 = (axis+1)%3;
    int axis2 = (axis+2)%3;

    // n0     grid size normal to face
    // n1,n2: grid size within face
    int n0 = index_global[axis0][1] - index_global[axis0][0];
    int n1 = index_global[axis1][1] - index_global[axis1][0];
    int n2 = index_global[axis2][1] - index_global[axis2][0];

    // index_global[][] should be divisible by r_factor_**level.  Just
    // test r_factor_ here.
    bool l0 = (index_global[axis1][0]/r_factor_)*r_factor_ == index_global[axis1][0];
    bool l1 = (index_global[axis1][1]/r_factor_)*r_factor_ == index_global[axis1][1];

    if (!l0) printf("index_global[%d][0] = %d\n",axis1,index_global[axis1][0]);
    assert(l0);
    if (!l1) printf("index_global[%d][1] = %d\n",axis1,index_global[axis1][1]);
    assert(l1);

    for (int face=0; face<2; face++) {

      // Loop over face zones that are aligned with coarse zones (hence "+= r")
      for (int index1=0; index1<n1; index1 += r_factor_) {
	for (int index2=0; index2<n2; index2 += r_factor_) {

	  AMRsolve_Grid* adjacent = grid.faces().adjacent(axis,face,index1,index2);

	  AMRsolve_Faces::Label& fz = grid.faces().label(axis,face,index1,index2);

	  // Add graph entries iff grid or adjacent grid is local, and if
	  // adjacent grid (if it exists) is in the next-coarser level
	  bool is_local =
	    ((adjacent != NULL) && (adjacent->is_local() || grid.is_local()));
	  bool is_coarse = (fz == AMRsolve_Faces::_coarse_);

	  if (is_local && is_coarse) {

	    // (fine) grid global indices
	    int index_fine[3]; 
	    index_fine[axis0] = index_global[axis0][0] + face*(n0 - r_factor_);
	    index_fine[axis1] = index_global[axis1][0] + index1;
	    index_fine[axis2] = index_global[axis2][0] + index2;

	    // (coarse) adjacent global indices
	    int index_coarse[3]; 
	    index_coarse[axis0] = (index_fine[axis0]) / r_factor_ + (face*r_factor_-1);
	    index_coarse[axis1] = (index_fine[axis1]) / r_factor_;
	    index_coarse[axis2] = (index_fine[axis2]) / r_factor_;

	    // adjust for periodicity
	    if (hierarchy_->is_periodic(axis0)) {
	      int period = hierarchy_->period_index(axis0,level_coarse);
 	      index_coarse[axis0] = (index_coarse[axis0] + period) % period;
 	    }

	    //--------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE 
	    //--------------------------------------------------

	    if (discret_type == discret_type_const) {

	      update_fine_coarse_const_(face,grid,axis0,phase,
					level_fine,level_coarse,
					index_fine,index_coarse);

	      //--------------------------------------------------
	      // GRAPH ENTRY: COARSE-TO-FINE
	      //--------------------------------------------------
	      
	      if (adjacent->is_local()) {

		update_coarse_fine_const_(face,*adjacent,axis0,phase,
					  level_fine,level_coarse,
					  index_fine,index_coarse);

	      }

	    } else {
	      char error_message[80];
	      strcpy(error_message,"Unknown parameter discret = ");
	      strcat(error_message,parameters_->value("discret").c_str());
	      ERROR(error_message);
	    } // if discret unexpected
	  } // if is_local && fz == AMRsolve_Faces::_coarse_
	} // for index2
      } // for index1
    } // for face
  } // for axis
} // AMRsolve_Hypre_FLD::init_nonstencil_()

//------------------------------------------------------------------------

/// Initialize matrix stencil and graph entries
void AMRsolve_Hypre_FLD::init_elements_matrix_()
{

  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level* level = itl++) {

    int part = level->index();

    // 1. Set stencil values within level
    ItLevelGridsLocal itlg (*level);
    while (AMRsolve_Grid* grid = itlg++)  init_matrix_stencil_(*grid);

    if (part > 0) {
      // *** WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all
      // *** grids; however, we only need to loop over "parent-child
      // *** pairs such that either child or parent is local to this
      // *** MPI process."
 
      // Set matrix values between levels
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid* grid = itag++)  init_matrix_nonstencil_(*grid);

    } // while level > 0
  } // while level = itl++

  // Clean up stencil connections between levels
  for (int part=1; part<hierarchy_->num_levels(); part++) 
    init_matrix_clear_(part);

} // init_elements_matrix_()

//------------------------------------------------------------------------

/// Set right-hand-side elements based on values in grids
void AMRsolve_Hypre_FLD::init_elements_rhs_()
{
  // declare shortcut variables
  Scalar Eavg, Ed_zl, Ed_yl, Ed_xl, Ed_xr, Ed_yr, Ed_zr, R, R0, kap, kap0;
  Scalar D_zl, D_yl, D_xl, D_xr, D_yr, D_zr;
  Scalar D0_zl, D0_yl, D0_xl, D0_xr, D0_yr, D0_zr;
  Scalar afac  = adot_  / aval_;
  Scalar afac0 = adot0_ / aval0_;
  Scalar dtfac  = dt_ * theta_;
  Scalar dtfac0 = dt_ * (1.0 - theta_);
  Scalar Rmin = 1.0e-20;
  Scalar c = 2.99792458e10;

  _TRACE_;
  // iterate over all grids local to this processor
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // for each local grid, use Enzo's "f" array to store RHS entries 
    // for the linear system.
    int n0,n1,n2;
    Scalar* values = grid->get_f(&n0,&n1,&n2);

    // access relevant arrays from this grid to compute RHS
    Scalar* E    = grid->get_E();
    Scalar* eta  = grid->get_eta();
    Scalar* HI   = grid->get_HI();
    Scalar* HeI  = grid->get_HeI();
    Scalar* HeII = grid->get_HeII();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_Ghosts(ghosts);
    int en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    int en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    int en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

    // access this grid's mesh spacing, and set relevant shortcuts
    Scalar hx,hy,hz;
    grid->h(hx,hy,hz);
    Scalar dxi  = aval_ / hx / lUn_;
    Scalar dyi  = aval_ / hy / lUn_;
    Scalar dzi  = aval_ / hz / lUn_;
    Scalar dxi0 = aval0_ / hx / lUn0_;
    Scalar dyi0 = aval0_ / hy / lUn0_;
    Scalar dzi0 = aval0_ / hz / lUn0_;

    // iterate over the domain, computing the RHS array
    int k2, k1, k0, k_l00, k_r00, k_0l0, k_0r0, k_00l, k_00r, k_000;
    int i0, i1, i2, i;
    for (i2=0; i2<n2; i2++) {
      k2 = ghosts[2][0] + i2;

      for (i1=0; i1<n1; i1++) {
	k1 = ghosts[1][0] + i1;

	for (i0=0; i0<n0; i0++) {
	  k0 = ghosts[0][0] + i0;

	  // compute indices of neighboring Enzo cells
	  k_l00 = k0-1 + en0*(k1   + en1*k2);
	  k_0l0 = k0   + en0*(k1-1 + en1*k2);
	  k_00l = k0   + en0*(k1   + en1*(k2-1));
	  k_000 = k0   + en0*(k1   + en1*k2);
	  k_r00 = k0+1 + en0*(k1   + en1*k2);
	  k_0r0 = k0   + en0*(k1+1 + en1*k2);
	  k_00r = k0   + en0*(k1   + en1*(k2+1));

	  //--------------
	  // z-directional limiter, lower face
	  Ed_zl = E[k_000] - E[k_00l];
	  Eavg = (E[k_000] + E[k_00l])*0.5;
	  R  = MAX(dzi *fabs(Ed_zl)/Eavg, Rmin);
	  R0 = MAX(dzi0*fabs(Ed_zl)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_00l])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_00l])*HIconst_ +
		   (HeI[k_000]  + HeI[k_00l])*HeIconst_ +
		   (HeII[k_000] + HeII[k_00l])*HeIIconst_)*0.5;
	  }
	  D_zl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_zl = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  //--------------
	  // y-directional limiter, lower face
	  Ed_yl = E[k_000] - E[k_0l0];
	  Eavg = (E[k_000] + E[k_0l0])*0.5;
	  R  = MAX(dyi *fabs(Ed_yl)/Eavg, Rmin);
	  R0 = MAX(dyi0*fabs(Ed_yl)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_0l0])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_0l0])*HIconst_ +
		   (HeI[k_000]  + HeI[k_0l0])*HeIconst_ +
		   (HeII[k_000] + HeII[k_0l0])*HeIIconst_)*0.5;
	  }
	  D_yl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_yl = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  //--------------
	  // x-directional limiter, lower face
	  Ed_xl = E[k_000] - E[k_l00];
	  Eavg = (E[k_000] + E[k_l00])*0.5;
	  R  = MAX(dxi *fabs(Ed_xl)/Eavg, Rmin);
	  R0 = MAX(dxi0*fabs(Ed_xl)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_l00])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_l00])*HIconst_ +
		   (HeI[k_000]  + HeI[k_l00])*HeIconst_ +
		   (HeII[k_000] + HeII[k_l00])*HeIIconst_)*0.5;
	  }
	  D_xl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_xl = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  //--------------
	  // x-directional limiter, upper face
	  Ed_xr = E[k_r00] - E[k_000];
	  Eavg = (E[k_r00] + E[k_000])*0.5;
	  R  = MAX(dxi *fabs(Ed_xr)/Eavg, Rmin);
	  R0 = MAX(dxi0*fabs(Ed_xr)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_r00])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_r00])*HIconst_ +
		   (HeI[k_000]  + HeI[k_r00])*HeIconst_ +
		   (HeII[k_000] + HeII[k_r00])*HeIIconst_)*0.5;
	  }
	  D_xr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_xr = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  //--------------
	  // y-directional limiter, upper face
	  Ed_yr = E[k_0r0] - E[k_000];
	  Eavg = (E[k_0r0] + E[k_000])*0.5;
	  R  = MAX(dyi *fabs(Ed_yr)/Eavg, Rmin);
	  R0 = MAX(dyi0*fabs(Ed_yr)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_0r0])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_0r0])*HIconst_ +
		   (HeI[k_000]  + HeI[k_0r0])*HeIconst_ +
		   (HeII[k_000] + HeII[k_0r0])*HeIIconst_)*0.5;
	  }
	  D_yr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_yr = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  //--------------
	  // z-directional limiter, upper face
	  Ed_zr = E[k_00r] - E[k_000];
	  Eavg = (E[k_00r] + E[k_000])*0.5;
	  R  = MAX(dzi *fabs(Ed_zr)/Eavg, Rmin);
	  R0 = MAX(dzi0*fabs(Ed_zr)/Eavg, Rmin);
	  if (Nchem_ == 1) {
	    kap = (HI[k_000] + HI[k_00r])*HIconst_*0.5;
	  } else {
	    kap = ((HI[k_000]   + HI[k_00r])*HIconst_ +
		   (HeI[k_000]  + HeI[k_00r])*HeIconst_ +
		   (HeII[k_000] + HeII[k_00r])*HeIIconst_)*0.5;
	  }
	  D_zr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	  D0_zr = c/sqrt(9.0*kap*kap*nUn0_*nUn0_ + R0*R0);

	  // opacity values in this cell
	  if (Nchem_ == 1) {
	    kap = HI[k_000]*HIconst_;
	  } else {
	    kap = HI[k_000]*HIconst_ + HeI[k_000]*HeIconst_ + HeII[k_000]*HeIIconst_;
	  }
	  kap0 = kap*nUn0_;
	  kap *= nUn_;

	  // set rhs in this cell
	  i = i0 + n0*(i1 + n1*i2);
	  values[i] = ( (dtfac/rUn_ + dtfac0/rUn0_)*eta[k_000]
		      + (1.0 - dtfac0*(afac0+c*kap0))*E[k_000]
		      + dtfac0*dxi0*dxi0*(D0_xr*Ed_xr-D0_xl*Ed_xl)
		      + dtfac0*dyi0*dyi0*(D0_yr*Ed_yr-D0_yl*Ed_yl)
		      + dtfac0*dzi0*dzi0*(D0_zr*Ed_zr-D0_zl*Ed_zl)
		      - (1.0 + dtfac*(afac+c*kap))*E[k_000]
		      + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)
		      + dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)
		      + dtfac*dzi*dzi*(D_zr*Ed_zr-D_zl*Ed_zl) );
	}
      }
    }

    // Set Hypre B_ vector to RHS values
    int part = grid->level();
    int lower[3] = { grid->index_lower(0), 
		     grid->index_lower(1), 
		     grid->index_lower(2) };
    int upper[3] = { grid->index_upper(0), 
		     grid->index_upper(1), 
		     grid->index_upper(2) };
    HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,values);

  }  // while grid = itg++

  //// WHY DO WE NOT CALL HYPRE_SStructFACZeroAMRVectorData HERE???? ////
  

} // init_elements_rhs_()

//------------------------------------------------------------------------

/// Compute a relative difference norm between E and E0, using the pnorm
/// and absolute tolerances given as input.  This routine ensures that 
/// for (pnorm != 0), the vector values are scaled by the cell volume, 
/// and that overlapped cells are not double-counted.
double AMRsolve_Hypre_FLD::rdiff_norm(double pnorm, double atol)
{

  // local variables
  int n0, n1, n2, en0, en1, en2, ghosts[3][2], part, lower[3], upper[3];
  int i, i0, i1, i2, k;
  Scalar *E, *E0, *values;
  Scalar hx, hy, hz, dV, weight, diff;
  double tmp, loc_diff, proc_diff=0.0, glob_diff=0.0;

  _TRACE_;
  // for each local grid, fill the "b" array to store difference values 
  //    iterate over all grids local to this processor
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // access AMRsolve array
    values = grid->get_f(&n0,&n1,&n2);

    // access relevant arrays from this grid to compute RHS
    E  = grid->get_E();
    E0 = grid->get_E0();

    // get buffering information on relating amrsolve grid to Enzo data
    grid->get_Ghosts(ghosts);
    en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

    // access this grid's mesh spacing, and set relevant shortcuts
    grid->h(hx,hy,hz);
    dV = (pnorm > 0.0) ? hx*hy*hz : 1.0;

    // iterate over the domain, computing the difference array
    for (i2=0; i2<n2; i2++)
      for (i1=0; i1<n1; i1++)
	for (i0=0; i0<n0; i0++) {

	  // compute Enzo cell index (k), AMRsolve grid index (i)
	  k = (ghosts[0][0] + i0) + en0*((ghosts[1][0] + i1) + en1*(ghosts[2][0] + i2));
	  i = i0 + n0*(i1 + n1*i2);

	  // compute local weight factor
	  weight = sqrt(fabs(E[k]*E0[k]) + atol);

	  // compute local weighted difference, and store inside b
	  values[i] = fabs(E[k] - E0[k])*dV/weight;
	}

    // Set Hypre B_ vector to RHS values
    part = grid->level();
    grid->get_limits(lower,upper);
    HYPRE_SStructVectorSetBoxValues(B_,part,lower,upper,0,values);
    
  }  // while grid = itg++

  // call HYPRE to zero out overlapped cell values
  int nlevels = hierarchy_->num_levels();
  if (nlevels > 1) {
    int plevels[nlevels];
    int rfactors[nlevels][3];
    for (int part=0; part<nlevels; part++) {
      plevels[part] = part;
      rfactors[part][0] = r_factor_;
      rfactors[part][1] = r_factor_;
      rfactors[part][2] = r_factor_;
    }
    HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
  }

  
  // iterate back over local grids, accumulating the difference norm
  ItHierarchyGridsLocal itg2(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // get level, grid information
    int level = grid->level();
    grid->get_limits(lower,upper);

    // extract Enzo solution
    values = grid->get_f(&n0,&n1,&n2);
    HYPRE_SStructVectorGather(B_);
    HYPRE_SStructVectorGetBoxValues(B_, level, lower, upper, 0, values);  

    // iterate over the domain, updating local p-norm 
    loc_diff = 0.0;
    if (pnorm > 0.0) { 
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  for (i0=0; i0<n0; i0++) {
	    tmp = values[i0 + n0*(i1 + n1*i2)];
	    loc_diff += pow(tmp,pnorm);
	  }
      proc_diff += loc_diff;
    // iterate over the domain, updating local max norm 
    } else { 
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  for (i0=0; i0<n0; i0++) {
	    tmp = values[i0 + n0*(i1 + n1*i2)];
	    loc_diff = (loc_diff > tmp) ? loc_diff : tmp;
	  }
      proc_diff = (loc_diff > proc_diff) ? loc_diff : proc_diff;
    }
   
  } // grid = itg2++

  // communicate to obtain overall sum/max
  if (pnorm > 0.0) {
    if (MPI_Allreduce(&proc_diff,&glob_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) != 0) {
      char error_message[80];
      strcpy(error_message,"rdiff_norm: Error in MPI_Allreduce!\n");
      ERROR(error_message);
    }
    glob_diff = pow(glob_diff, 1.0/pnorm);
  } else {
    if (MPI_Allreduce(&proc_diff,&glob_diff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) != 0) {
      char error_message[80];
      strcpy(error_message,"rdiff_norm: Error in MPI_Allreduce!\n");
      ERROR(error_message);
    }
  }

  // return with overall difference value
  return glob_diff;

} // rdiff_norm()

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid interior
void AMRsolve_Hypre_FLD::init_matrix_stencil_(AMRsolve_Grid& grid)
{
  int n          = grid.num_unknowns();
  int entries[7] = { 0,1,2,3,4,5,6 };
  int    n3[3]   = {grid.n(0),grid.n(1),grid.n(2)};

  double* v0;         // Diagonal elements
  double* v1[3][2];   // Off-diagonal elements

  // declare shortcut variables
  double Eavg, Ed_zl, Ed_yl, Ed_xl, Ed_xr, Ed_yr, Ed_zr, R, kap;
  double D_zl, D_yl, D_xl, D_xr, D_yr, D_zr;
  double afac = adot_ / aval_;
  double dtfac = dt_ * theta_;
  double Rmin = 1.0e-20;
  double c = 2.99792458e10;

  // access relevant arrays from this grid to compute RHS
  Scalar* E    = grid.get_E();
  Scalar* HI   = grid.get_HI();
  Scalar* HeI  = grid.get_HeI();
  Scalar* HeII = grid.get_HeII();
  
  // get buffering information on relating amrsolve grid to Enzo data
  int ghosts[3][2]; 
  grid.get_Ghosts(ghosts);
  int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
  int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
  int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

  // access this grid's mesh spacing, and set relevant shortcuts
  double dxi = aval_ / grid.h(0) / lUn_;
  double dyi = aval_ / grid.h(1) / lUn_;
  double dzi = aval_ / grid.h(2) / lUn_;
  
  // Allocate storage
  v0 = new double[n];
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      v1[axis][face] = new double[n];
    } // for face
  } // for axis

  //-----------------------------------------------------------
  // Set stencil for all unknowns, ignoring boundary conditions
  //-----------------------------------------------------------

  int k2, k1, k0, k_l00, k_r00, k_0l0, k_0r0, k_00l, k_00r, k_000;
  int i0, i1, i2, i;
  for (i2=0; i2<n3[2]; i2++) {
    k2 = ghosts[2][0] + i2;

    for (i1=0; i1<n3[1]; i1++) {
      k1 = ghosts[1][0] + i1;

      for (i0=0; i0<n3[0]; i0++) {
	k0 = ghosts[0][0] + i0;
	  
	// compute indices of neighboring Enzo cells
	k_l00 = k0-1 + en0*(k1   + en1*k2);
	k_0l0 = k0   + en0*(k1-1 + en1*k2);
	k_00l = k0   + en0*(k1   + en1*(k2-1));
	k_000 = k0   + en0*(k1   + en1*k2);
	k_r00 = k0+1 + en0*(k1   + en1*k2);
	k_0r0 = k0   + en0*(k1+1 + en1*k2);
	k_00r = k0   + en0*(k1   + en1*(k2+1));

	//--------------
	// z-directional limiter, lower face
	Ed_zl = E[k_000] - E[k_00l];
	Eavg = (E[k_000] + E[k_00l])*0.5;
	R = MAX(dzi*fabs(Ed_zl)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_00l])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_00l])*HIconst_ +
		 (HeI[k_000]  + HeI[k_00l])*HeIconst_ +
		 (HeII[k_000] + HeII[k_00l])*HeIIconst_)*0.5;
	}
	D_zl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	//--------------
	// y-directional limiter, lower face
	Ed_yl = E[k_000] - E[k_0l0];
	Eavg = (E[k_000] + E[k_0l0])*0.5;
	R = MAX(dyi*fabs(Ed_yl)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_0l0])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_0l0])*HIconst_ +
		 (HeI[k_000]  + HeI[k_0l0])*HeIconst_ +
		 (HeII[k_000] + HeII[k_0l0])*HeIIconst_)*0.5;
	}
	D_yl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	//--------------
	// x-directional limiter, lower face
	Ed_xl = E[k_000] - E[k_l00];
	Eavg = (E[k_000] + E[k_l00])*0.5;
	R = MAX(dxi*fabs(Ed_xl)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_l00])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_l00])*HIconst_ +
		 (HeI[k_000]  + HeI[k_l00])*HeIconst_ +
		 (HeII[k_000] + HeII[k_l00])*HeIIconst_)*0.5;
	}
	D_xl = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	//--------------
	// x-directional limiter, upper face
	Ed_xr = E[k_r00] - E[k_000];
	Eavg = (E[k_r00] + E[k_000])*0.5;
	R = MAX(dxi*fabs(Ed_xr)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_r00])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_r00])*HIconst_ +
		 (HeI[k_000]  + HeI[k_r00])*HeIconst_ +
		 (HeII[k_000] + HeII[k_r00])*HeIIconst_)*0.5;
	}
	D_xr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	//--------------
	// y-directional limiter, upper face
	Ed_yr = E[k_0r0] - E[k_000];
	Eavg = (E[k_0r0] + E[k_000])*0.5;
	R = MAX(dyi*fabs(Ed_yr)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_0r0])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_0r0])*HIconst_ +
		 (HeI[k_000]  + HeI[k_0r0])*HeIconst_ +
		 (HeII[k_000] + HeII[k_0r0])*HeIIconst_)*0.5;
	}
	D_yr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	//--------------
	// z-directional limiter, upper face
	Ed_zr = E[k_00r] - E[k_000];
	Eavg = (E[k_00r] + E[k_000])*0.5;
	R = MAX(dzi*fabs(Ed_zr)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_000] + HI[k_00r])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_000]   + HI[k_00r])*HIconst_ +
		 (HeI[k_000]  + HeI[k_00r])*HeIconst_ +
		 (HeII[k_000] + HeII[k_00r])*HeIIconst_)*0.5;
	}
	D_zr = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);
	
	// opacity values in this cell
	if (Nchem_ == 1) {
	  kap = HI[k_000]*HIconst_;
	} else {
	  kap = HI[k_000]*HIconst_ + HeI[k_000]*HeIconst_ + HeII[k_000]*HeIIconst_;
	}
	kap *= nUn_;
	
	// get the linear grid index
	i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);

	// set the matrix entries
	v1[2][0][i] = -dtfac*dzi*dzi*D_zl;    // z-left
	v1[1][0][i] = -dtfac*dyi*dyi*D_yl;    // y-left
	v1[0][0][i] = -dtfac*dxi*dxi*D_xl;    // x-left
	v1[0][1][i] = -dtfac*dxi*dxi*D_xr;    // x-right
	v1[1][1][i] = -dtfac*dyi*dyi*D_yr;    // y-right
	v1[2][1][i] = -dtfac*dzi*dzi*D_zr;    // z-right
	v0[i] = 1.0 + dtfac*(afac + c*kap + dxi*dxi*(D_xl+D_xr) 
		    + dyi*dyi*(D_yl+D_yr) + dzi*dzi*(D_zl+D_zr));  // self

      } // for i0
    } // for i1
  } // for i2

//   WHERE; printf("v0[0]=%g\n",v0[0]);

  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries
  //-----------------------------------------------------------

  // update matrix/rhs based on boundary conditions/location
  //    z-left face
  if (grid.faces().label(2,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[2][0] == 1) {           // Dirichlet
      i2 = 0;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[2][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[2][0] == 2) {    // Neumann
      i2 = 0;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][0][i];
	  v1[2][0][i] = 0.0;
	} 
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    y-left face
  if (grid.faces().label(1,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[1][0] == 1) {           // Dirichlet
      i1 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[1][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[1][0] == 2) {    // Neumann
      i1 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][0][i];
	  v1[1][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    x-left face
  if (grid.faces().label(0,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[0][0] == 1) {           // Dirichlet
      i0 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[0][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[0][0] == 2) {    // Neumann
      i0 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][0][i];
	  v1[0][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    x-right face
  if (grid.faces().label(0,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[0][1] == 1) {           // Dirichlet
      i0 = n3[0]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[0][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[0][1] == 2) {    // Neumann
      i0 = n3[0]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][1][i];
	  v1[0][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    y-right face
  if (grid.faces().label(1,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[1][1] == 1) {           // Dirichlet
      i1 = n3[1]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[1][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[1][1] == 2) {    // Neumann
      i1 = n3[1]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][1][i];
	  v1[1][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    z-right face
  if (grid.faces().label(2,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[2][1] == 1) {           // Dirichlet
      i2 = n3[2]-1;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[2][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[2][1] == 2) {    // Neumann
      i2 = n3[2]-1;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][1][i];
	  v1[2][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary
  

  //-----------------------------------------------------------
  // insert matrix entries into Hypre matrix A_
  //-----------------------------------------------------------

  //   AMRsolve_Faces& faces = grid.faces();
  int level = grid.level();
  int index_lower[3] = { grid.index_lower(0), 
			 grid.index_lower(1), 
			 grid.index_lower(2) };
  int index_upper[3] = { grid.index_upper(0), 
			 grid.index_upper(1), 
			 grid.index_upper(2) };

  // Update matrix stencil values
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper,
				  0, 1, &entries[0], v0);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[1], v1[0][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[2], v1[0][0]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[3], v1[1][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[4], v1[1][0]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[5], v1[2][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[6], v1[2][0]);

  // Deallocate arrays
  delete [] v0;
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] v1[axis][face];
    } // for face
  } // for axis

} // AMRsolve_Hypre_FLD::init_matrix_stencil_()

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver
void AMRsolve_Hypre_FLD::init_matrix_clear_(int part)
{
  if (part > 0) {
    int r_factors[3] = {r_factor_,r_factor_,r_factor_}; 
    HYPRE_SStructFACZeroAMRMatrixData(A_, part-1, r_factors);
  }

} // AMRsolve_Hypre_FLD::init_matrix_clear_()

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver
void AMRsolve_Hypre_FLD::solve_pfmg_(int itmax, double restol)
{
  _TRACE_;
  // Create and initialize the solver
  HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string srlxtype = parameters_->value("solver_rlxtype");
  std::string snpre    = parameters_->value("solver_npre");
  std::string snpost   = parameters_->value("solver_npost");
  std::string sprintl  = parameters_->value("solver_printl");
  std::string slog     = parameters_->value("solver_log");

  //   if not defined, then define them
  if (srlxtype == "")  srlxtype = "1";
  if (snpre == "")     snpre = "1";
  if (snpost == "")    snpost = "1";
  if (sprintl == "")   sprintl = "1";
  if (slog == "")      slog = "1";

  //   set local variables
  int rlxtype = atoi(srlxtype.c_str());
  int npre    = atoi(snpre.c_str());
  int npost   = atoi(snpost.c_str());
  int printl  = atoi(sprintl.c_str());
  int log     = atoi(slog.c_str());

  // set solver options
  if (itmax != 0 )   HYPRE_SStructSysPFMGSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructSysPFMGSetTol(solver_,    restol);
  HYPRE_SStructSysPFMGSetRelaxType(solver_,    rlxtype);
  HYPRE_SStructSysPFMGSetNumPreRelax(solver_,  npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(solver_, npost);
  HYPRE_SStructSysPFMGSetPrintLevel(solver_,   printl);
  HYPRE_SStructSysPFMGSetLogging(solver_,      log);

  // setup solver 
  HYPRE_SStructSysPFMGSetup(solver_,A_,B_,X_);

  // Solve the linear system
  HYPRE_SStructSysPFMGSolve(solver_,A_,B_,X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructSysPFMGGetNumIterations(solver_,&iter_);
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver_,&resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre PFMG num iterations: %d\n",iter_);
    printf("hypre PFMG final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructSysPFMGDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_FLD::solve_pfmg_()

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver
void AMRsolve_Hypre_FLD::solve_fac_(int itmax, double restol)
{
  _TRACE_;
  int i;

  // Create the solver
  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);

  // Initialize parts
  int num_parts = hierarchy_->num_levels();
  HYPRE_SStructFACSetMaxLevels(solver_, num_parts);
  int *parts = new int[num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;
  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);

  // Initialize refinement factors
  int3 *refinements = new int3[num_parts];
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r_factor_;
    refinements[i][1] = r_factor_;
    refinements[i][2] = r_factor_;
  }
  HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);

  // extract some additional solver parameters
  std::string srlxtype = parameters_->value("solver_rlxtype");
  std::string snpre    = parameters_->value("solver_npre");
  std::string snpost   = parameters_->value("solver_npost");
  std::string scsolve  = parameters_->value("solver_csolve");
  std::string sprintl  = parameters_->value("solver_printl");
  std::string slog     = parameters_->value("solver_log");

  //   if not defined, then define them
  if (srlxtype == "")  srlxtype = "2";
  if (snpre == "")     snpre = "2";
  if (snpost == "")    snpost = "2";
  if (scsolve == "")   scsolve = "1";
  if (sprintl == "")   sprintl = "1";
  if (slog == "")      slog = "1";

  //   set local variables
  int relax   = atoi(srlxtype.c_str());
  int npre    = atoi(snpre.c_str());
  int npost   = atoi(snpost.c_str());
  int printl  = atoi(sprintl.c_str());
  int log     = atoi(slog.c_str());
  int csolve  = atoi(scsolve.c_str());

  // solver parameters
  HYPRE_SStructFACSetNumPreRelax(solver_,      npre);
  HYPRE_SStructFACSetNumPostRelax(solver_,     npost);
  HYPRE_SStructFACSetCoarseSolverType(solver_, csolve);
  HYPRE_SStructFACSetRelaxType(solver_,        relax);

  // stopping criteria
  if (itmax != 0 )    HYPRE_SStructFACSetMaxIter(solver_, itmax);
  if (restol != 0.0)  HYPRE_SStructFACSetTol(solver_,     restol);

  // output amount
  HYPRE_SStructFACSetLogging(solver_, 1);

  // prepare for solve
  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructFACGetNumIterations(solver_, &iter_);
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre FAC num iterations: %d\n",iter_);
    printf("hypre FAC final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructFACDestroy2(solver_);

  // Delete local dynamic storage
  delete [] parts;
  delete [] refinements;

  _TRACE_;
} // AMRsolve_Hypre_FLD::solve_fac_()

//------------------------------------------------------------------------

/// Initialize the BICGSTAB hypre solver
void AMRsolve_Hypre_FLD::solve_bicgstab_(int itmax, double restol)
{
  _TRACE_;
  // Create the solver
  HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 )   HYPRE_SStructBiCGSTABSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructBiCGSTABSetTol(solver_,    restol);

  // output amount
  HYPRE_SStructBiCGSTABSetLogging(solver_, log);

  // Initialize the solver
  HYPRE_SStructBiCGSTABSetup(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructBiCGSTABSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
  HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructBiCGSTABDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_FLD::solve_bicgstab_()

//------------------------------------------------------------------------

/// Initialize the BiCGSTAB hypre solver with BoomerAMG preconditioning
void AMRsolve_Hypre_FLD::solve_bicgstab_boomer_(int itmax, double restol)
{
  _TRACE_;
  HYPRE_Solver solver;

  // Create the solver
  HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);

  // stopping criteria
  if (itmax != 0 )   HYPRE_BiCGSTABSetMaxIter(solver,itmax);
  if (restol != 0.0) HYPRE_BiCGSTABSetTol(solver,    restol);

  // extract some additional solver parameters
  std::string slog    = parameters_->value("solver_log");
  std::string sprintl = parameters_->value("solver_printl");
  if (slog == "")      slog = "1";
  if (sprintl == "")   sprintl = "1";
  int log    = atoi(slog.c_str());
  int printl = atoi(sprintl.c_str());

  // output amount
  HYPRE_BiCGSTABSetLogging(solver, log);

  // Set BoomerAMG preconditioner
  HYPRE_Solver par_precond;
  HYPRE_BoomerAMGCreate(&par_precond);
  HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
  HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
  HYPRE_BoomerAMGSetTol(par_precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(par_precond, printl);
  HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
  HYPRE_BoomerAMGSetMaxIter(par_precond, 1);
  HYPRE_BiCGSTABSetPrecond(solver,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
			   par_precond);

  // Initialize the solver
  HYPRE_BiCGSTABSetup(solver, (HYPRE_Matrix) A_, 
		      (HYPRE_Vector) B_, (HYPRE_Vector) X_);

  // Solve the linear system
  HYPRE_BiCGSTABSolve(solver, (HYPRE_Matrix) A_, 
		      (HYPRE_Vector)B_, (HYPRE_Vector) X_);

  // Write out some diagnostic info about the solve
  HYPRE_BiCGSTABGetNumIterations(solver, &iter_);
  HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &resid_);
  if (debug &&pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }


  // Delete the solver and preconditiner
  HYPRE_BoomerAMGDestroy(par_precond);
  HYPRE_BiCGSTABDestroy(solver);

  _TRACE_;
} // AMRsolve_Hypre_FLD::solve_bicgstab_boomer_()

//------------------------------------------------------------------------

/// Initialize the GMRES hypre solver
void AMRsolve_Hypre_FLD::solve_gmres_(int itmax, double restol)
{
  _TRACE_;
  // Create the solver
  HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 )   HYPRE_SStructGMRESSetMaxIter(solver_, itmax);
  if (restol != 0.0) HYPRE_SStructGMRESSetTol(solver_,     restol);

  // output amount
  HYPRE_SStructGMRESSetLogging(solver_, log);

  // Initialize the solver
  HYPRE_SStructGMRESSetup(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructGMRESSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructGMRESGetNumIterations(solver_, &iter_);
  HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre GMRES num iterations: %d\n",iter_);
    printf("hypre GMRES final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructGMRESDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_FLD::solve_gmres_()

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the 
/// fine-grid matrix to account for the coarse grid neighbor.  Called in 
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_Hypre_FLD::update_fine_coarse_const_(int face, 
						   AMRsolve_Grid& grid_fine, 
						   int axis0, 
						   phase_enum phase,
						   int level_fine, 
						   int level_coarse,
						   int index_fine[3], 
						   int index_coarse[3])
{
  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  // declare shortcut variables
  Scalar Eavg, Ed, R, kap, D;
  Scalar afac = adot_ / aval_;
  Scalar dtfac = dt_ * theta_;
  Scalar Rmin = 1.0e-20;
  Scalar c = 2.99792458e10;

  // access relevant arrays from this grid to compute RHS
  Scalar* E    = grid_fine.get_E();
  Scalar* HI   = grid_fine.get_HI();
  Scalar* HeI  = grid_fine.get_HeI();
  Scalar* HeII = grid_fine.get_HeII();
  
  // get active enzo grid size
  int n3[3];
  grid_fine.get_size(n3);

  // get buffering information on relating amrsolve grid to Enzo data
  int ghosts[3][2]; 
  grid_fine.get_Ghosts(ghosts);
  int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
  int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
  int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

  // indices for Enzo fine grid cells on both sides of interface
  int k_row, k_col, k0, k1, k2;
  int adj0=0, adj1=0, adj2=0;
  if (axis0 == 0)  adj0 = (face == 0) ? -1 : 1;  
  if (axis0 == 1)  adj1 = (face == 0) ? -1 : 1;  
  if (axis0 == 2)  adj2 = (face == 0) ? -1 : 1;  
  
  // set this grid's mesh spacing in this direction
  Scalar dxi;
  if (axis0 == 0)  dxi = aval_ / grid_fine.h(0) / lUn_;
  if (axis0 == 1)  dxi = aval_ / grid_fine.h(1) / lUn_;
  if (axis0 == 2)  dxi = aval_ / grid_fine.h(2) / lUn_;
  
  //--------------------------------------------------
  // (*) CONSTANT
  //     Scale        = 2/3
  //     Coefficients = determined on the fly
  //--------------------------------------------------

  int index_increment[][3] = {{face*(r_factor_-1),0,0},
			      {0,1,0},
			      {0,0,1},
			      {0,-1,0},
			      {-face*(r_factor_-1),0,-1}};

  if (grid_fine.is_local()) {

    if (phase == phase_graph) {

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {
	HYPRE_SStructGraphAddEntries(graph_, level_fine, index_fine, 
				     0, level_coarse, index_coarse, 0);
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];
      } // for k = 1:4

    } else if (phase == phase_matrix) {

      // fine->coarse off-diagonal scaling
      double val_s = 2.0 / 3.0;
      int entry;
      double val, value;

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {

	// set indices for Enzo fine cells on both sides of interface
	k0 = index_fine[0] + ghosts[0][0];
	k1 = index_fine[1] + ghosts[1][0];
	k2 = index_fine[2] + ghosts[2][0];
	k_row = k0 + en0*(k1 + en1*k2);
	k_col = k0 + adj0 + en0*(k1 + adj1 + en1*(k2 + adj2));

	// Compute limiter at this face
	Ed = E[k_row] - E[k_col];
	Eavg = (E[k_row] + E[k_col])*0.5;
	R = MAX(dxi*fabs(Ed)/Eavg, Rmin);
	if (Nchem_ == 1) {
	  kap = (HI[k_row] + HI[k_col])*HIconst_*0.5;
	} else {
	  kap = ((HI[k_row]   + HI[k_col])*HIconst_ +
		 (HeI[k_row]  + HeI[k_col])*HeIconst_ +
		 (HeII[k_row] + HeII[k_col])*HeIIconst_)*0.5;
	}
	D = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);


	// Set matrix values across coarse/fine face
	val = -val_s * dtfac * dxi * dxi * D;
	
	//   Update off-diagonal
	entry = grid_fine.counter(index_fine)++;
	value = val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	//   Update diagonal
	entry = 0;
	value = -val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);


	// Clear original matrix values from stencil
	val = -dtfac * dxi * dxi * D;
	
	//   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
	entry = 2*axis0 + 1 + (1-face);
	value = -val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	//   Update diagonal
	entry = 0;
	value = val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	// Update indices
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];

      } // for k = 1,4
    } // if phase == phase_matrix
  } // if grid_fine.is_local()
} // if discret_type_const

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the 
/// coarse-grid matrix to account for the fine grid neighbors.  Called in 
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_Hypre_FLD::update_coarse_fine_const_(int face, 
						   AMRsolve_Grid& grid_coarse, 
						   int axis0, 
						   phase_enum phase,
						   int level_fine, 
						   int level_coarse,
						   int index_fine[3], 
						   int index_coarse[3])
{
  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  // set graph entry
  if (phase == phase_graph) {

    int index_increment[][3] = {{1,0,0},
				{0,1,0},
				{-1,0,0},
				{0,0,1},
				{1,0,0},
				{0,-1,0},
				{-1,0,0},
				{0,0,-1}};

    for (int k=0; k<8; k++) {
      HYPRE_SStructGraphAddEntries(graph_, level_coarse, index_coarse, 
				   0, level_fine, index_fine, 0);
      index_fine[0] += index_increment[k][0];
      index_fine[1] += index_increment[k][1];
      index_fine[2] += index_increment[k][2];
    } // for k=0,7

  // set matrix entry
  } else if (phase == phase_matrix) {

    // declare shortcut variables
    double Eavg, Ed, R, kap, D;
    double afac = adot_ / aval_;
    double dtfac = dt_ * theta_;
    double Rmin = 1.0e-20;
    double c = 2.99792458e10;
    
    // access relevant arrays from this grid to compute RHS
    Scalar* E    = grid_coarse.get_E();
    Scalar* HI   = grid_coarse.get_HI();
    Scalar* HeI  = grid_coarse.get_HeI();
    Scalar* HeII = grid_coarse.get_HeII();
    
    // get active enzo grid size
    int n3[3];
    grid_coarse.get_size(n3);

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid_coarse.get_Ghosts(ghosts);
    int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions
    
    // set this grid's mesh spacing in this direction
    Scalar dxi;
    if (axis0 == 0)  dxi = aval_ / grid_coarse.h(0) / lUn_;
    if (axis0 == 1)  dxi = aval_ / grid_coarse.h(1) / lUn_;
    if (axis0 == 2)  dxi = aval_ / grid_coarse.h(2) / lUn_;
    
    // set indices for Enzo coarse grid cells on both sides of interface
    int adj0=0, adj1=0, adj2=0;
    if (axis0 == 0)  adj0 = (face == 0) ? -1 : 1;  
    if (axis0 == 1)  adj1 = (face == 0) ? -1 : 1;  
    if (axis0 == 2)  adj2 = (face == 0) ? -1 : 1;  
    
    int k0 = index_coarse[0] + ghosts[0][0];
    int k1 = index_coarse[1] + ghosts[1][0];
    int k2 = index_coarse[2] + ghosts[2][0];
    int k_row = k0 + en0*(k1 + en1*k2);
    int k_col = k0 + adj0 + en0*(k1 + adj1 + en1*(k2 + adj2));
    
    // Compute limiter at this face
    Ed = E[k_row] - E[k_col];
    Eavg = (E[k_row] + E[k_col])*0.5;
    R = MAX(dxi*fabs(Ed)/Eavg, Rmin);
    if (Nchem_ == 1) {
      kap = (HI[k_row] + HI[k_col])*HIconst_*0.5;
    } else {
      kap = ((HI[k_row]   + HI[k_col])*HIconst_ +
	     (HeI[k_row]  + HeI[k_col])*HeIconst_ +
	     (HeII[k_row] + HeII[k_col])*HeIIconst_)*0.5;
    }
    D = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);

    // set matrix entry across coarse/fine face
    double val_s = 1./8.;
    double val = -val_s * dtfac * dxi * dxi * D;
    int    entry;
    double value;

    // Adjust coarse-fine nonstencil values
    for (int i=0; i<8; i++) {
      // Set new nonstencil coarse-fine entry
      entry = grid_coarse.counter(index_coarse)++;
      value = val;
      HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				     0, 1, &entry, &value);
      // Adjust stencil diagonal
      entry = 0;
      value = -val;
      HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				     0, 1, &entry, &value);
    } // for i=0,7

    // Clear original matrix values from stencil
    val = -dtfac * dxi * dxi * D;

    //   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
    //   (note: "face" is for fine grid, but we want coarse)
    entry = 2*axis0 + 1 + face;
    value = -val;
    HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				   0, 1, &entry, &value);

    //   Update diagonal
    entry = 0;
    value = val;
    HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				   0, 1, &entry, &value);

  } // if phase == phase_matrix
}

//------------------------------------------------------------------------
