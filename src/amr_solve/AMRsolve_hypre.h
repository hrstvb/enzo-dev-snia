/// @file      newgrav-hypre.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Hypre class

#ifndef AMRSOLVE_HYPRE_H
#define AMRSOLVE_HYPRE_H

class AMRsolve_Hypre {

  /// @class    AMRsolve_Hypre
  /// @brief    AMRsolve_Hypre class for interfacing to the LLNL hypre solver

private:

  HYPRE_SStructGrid    grid_;          // hypre grid
  HYPRE_SStructGraph   graph_;         // hypre graph
  HYPRE_SStructStencil stencil_;       // hypre stencil
  HYPRE_SStructMatrix  A_;             // hypre matrix
  HYPRE_SStructVector  B_;             // hypre vector right-hand side
  HYPRE_SStructVector  X_;             // hypre vector solution
  HYPRE_SStructSolver  solver_;        // hypre solver

  AMRsolve_Parameters *parameters_;    // Pointer to parameters
  AMRsolve_Hierarchy  *hierarchy_;     // Pointer to the hierarchy

  double               resid_;         // Solver residual
  int                  iter_;          // Solver iterations

  const int            r_factor_;      // Refinement factor
  Scalar               matrix_scale_;  // 1.0:  1 1 1 -6 1 1 1

public:

  AMRsolve_Hypre(AMRsolve_Hierarchy& hierarchy, AMRsolve_Parameters& parameters);

  ~AMRsolve_Hypre();

  void init_hierarchy(AMRsolve_Mpi& mpi);
  void init_stencil();
  void init_graph();
  void init_elements(std::vector<AMRsolve_Point *> points, 
		     Scalar f_scale=1.0);
  void solve();
  void evaluate();

  int    iterations() { return iter_; };
  double residual() { return resid_; };

private:

  enum phase_enum  {phase_unknown,phase_graph,phase_matrix};

  // init_graph() functions
  void init_graph_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_graph); };

  // init_elements() functions
  void init_elements_matrix_();
  void init_elements_rhs_(std::vector<AMRsolve_Point *>& points, 
			  Scalar f_scale=1.0);

  void init_matrix_stencil_(AMRsolve_Grid& grid);
  void init_matrix_clear_(int part);
  void init_matrix_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_matrix); };

  void init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase);
  
  // init_vector() functions
  Scalar init_vector_points_(std::vector<AMRsolve_Point *>& points);
  Scalar init_vector_file_(std::string file_prefix, bool is_packed);
  Scalar init_vector_attach_(Scalar f_scale=1.0);

  // solve() functions
  void solve_fac_(int itmax, double restol);
  void solve_bicgstab_(int itmax, double restol);
  void solve_bicgstab_boomer_(int itmax, double restol);
  void solve_gmres_(int itmax, double restol);
  void solve_pfmg_(int itmax, double restol);

  // matrix graph update functions
  void update_fine_coarse_const_(int face, AMRsolve_Grid& grid, 
				 int axis0, phase_enum phase,
				 int level_fine, int level_coarse,
				 int igg3[3], int ign3[3]);

  void update_coarse_fine_const_(int face, AMRsolve_Grid& grid, 
				 int axis0, phase_enum phase,
				 int level_fine, int level_coarse,
				 int igg3[3], int ign3[3]);
};

#endif /* AMRSOLVE_HYPRE_H */
