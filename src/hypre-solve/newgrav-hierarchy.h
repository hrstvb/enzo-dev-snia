// $Id: newgrav-hierarchy.h 19 2010-03-19 23:04:28Z jbordner $

/// @file      newgrav-hierarchy.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Declaration of the Hierarchy class

#ifndef NEWGRAV_HIERRACHY_H
#define NEWGRAV_HIERRACHY_H

#ifdef HYPRE_GRAV
// ENZO DEPENDENCY
struct LevelHierarchyEntry;
#endif

#define LEVEL_UNKNOWN -1

class Hierarchy
{

  /// @class    Hierarchy
  /// @brief    Hierarchy class for representing an AMR hierarchy

  friend class ItHierarchyGridsLocal;
  friend class ItHierarchyGridsAll;
  friend class ItHierarchyLevels;
  friend class ItHierarchyLevelsReverse;

private:

  /// Dimension
  int                        dimension_;
  /// List of each grid's parent
  std::map<Grid *, Grid * >  grid_parent_;
  /// Global index of lowest vertex in globalo root grid
  int                 il0_[3];
  /// Number of zones in global root grid
  int                 n0_[3];

  /// Position of lowest vertex
  Scalar              xl_[3];
  /// Position of highest vertex
  Scalar              xu_[3];
  /// Domain periodicity, or 0 if not periodic
  double              period_domain_[3];
  /// Grid periodicity, or 0 if not periodic
  int                 period_index_[3];

protected:

  /// List of all grids, with extra sentinal 0-pointer at the end
  std::vector <Grid *>       grids0_;

  /// List of Levels,  with extra sentinal 0-pointer at the end
  std::vector <Level *>      levels0_;


 public:

  Hierarchy () throw ();

  ~Hierarchy () throw ();

#ifdef HYPRE_GRAV

  /// ENZO INTERFACE: attach to the Enzo hierarchy

  void enzo_attach (LevelHierarchyEntry *LevelArray[],
		    int level_coarse = LEVEL_UNKNOWN,
		    int level_fine   = LEVEL_UNKNOWN) throw ();

  /// ENZO INTERFACE: detach from the Enzo hierarchy

  void enzo_detach () throw ();

#endif

  /// Insert a Grid into the Hierarchy
  void insert_grid (Grid * grid) throw ();

  /// Initialize the Hierarchy given a Domain, Mpi, and whether it's periodic
  void initialize (Domain & domain, Mpi & mpi, bool is_periodic) throw();

  void set_dim (int d) throw () { dimension_ = d; };

  /// Write the hierarchy Levels and Grids to stdout
  void print () throw ();

  /// Write the Hierarchy to the given file.
  void write (FILE * fp = stdout) throw ();

  /// Write the Hierarchy grids to in geomview format
  void geomview_grids (Mpi & mpi) throw ();

  /// Checking validity of the Hierarchy [not implemented]

  void check () throw () {  
    printf ("Hierarchy::check() is not implemented yet\n");
  }

  /// Set a Grids parent
  void    set_parent (Grid * grid, Grid * parent) { grid_parent_[grid]=parent; };

  /// Return the dimension of the Hierarchy
  int dimension () throw () { return dimension_; } ;

  /// Return the ith Level of the Hierarchy.   No error checking.
  Level & level      (int i)            { return * levels0_.at(i); };

  /// Return the number of levels in the Hierarchy
  int     num_levels ()                 { return   levels0_.size() - 1; };

  /// Return the jth Grid.   No error checking.
  Grid &  return_grid       (int j)            { return * grids0_.at(j); };

  /// Return the number of grids
  int     num_grids  ()                 { return   grids0_.size() - 1; };

  /// Return the jth Grid in the ith Level of the Hierarchy.   No error checking.
  Grid &  return_grid       (int i, int j)     { return   levels0_[i]->return_grid(j); };

  /// Return the number of grids in the ith level.   No error checking.
  int     num_grids  (int i)            { return   levels0_[i]->num_grids(); };

  /// Return the Grids parent Grid
  Grid *  parent     (Grid & grid)      { return   grid_parent_[&grid]; };

  /// Return the lower index of the global root-grid.  No error checking on i.
  int i_lower0(int i) throw ()
  { return il0_[i]; };

  /// Return the lower indices of the global root-grid.
  void i_lower0(int &il0, int &il1, int &il2) throw ()
  { 
    il0 = il0_[0];
    il1 = il0_[1];
    il2 = il0_[2];
  };

  /// Return the upper index of the global root-grid.  No error checking on i.
  int i_upper0(int i) throw ()
  { return il0_[i] + n0_[i] - 1; };

  /// Return the upper indices of the global root-grid.
  void i_upper0(int &iu0, int &iu1, int &iu2) throw ()
  { 
    iu0 = il0_[0] + n0_[0] - 1;
    iu1 = il0_[1] + n0_[1] - 1;
    iu2 = il0_[2] + n0_[2] - 1;
  };

  /// Return lower and upper global indices of the root-grid.  Upper
  /// indices are incremented by one.

  void indices0(int ind[3][2]) throw ()
  { 
    ind[0][0] = il0_[0];
    ind[1][0] = il0_[1];
    ind[2][0] = il0_[2];
    ind[0][1] = il0_[0] + n0_[0];
    ind[1][1] = il0_[1] + n0_[1];
    ind[2][1] = il0_[2] + n0_[2];
  };


  /// Return the number of unknowns along the ith coordinate of the root grid.  No error checking on i.
  int num_unknowns0(int i) throw ()
  { return n0_[i]; };


  /// Return the total number of unknowns in the root grid.
  int num_unknowns0() throw ()
  { return n0_[0]*n0_[1]*n0_[2]; };

  /// Return periodicity of indices of the given level.
  int period_index(int axis, int level) throw ()
  { 
    int p=period_index_[axis];
    for (int i=0; i<level; i++) {
      p *= const_r_factor;
    }
    return p;
  }

  /// Return periodicity of indices of the given level.
  void period_index(int period3[3], int level) throw ()
  { 
    for (int i=0; i<3; i++) {
      period3[i] = period_index_[i];
      for (int l=0; l<level; l++) {
	period3[i] *= const_r_factor;
      }
    }
    return;
  }

  /// Return whether axis is periodic or not
  bool is_periodic(int axis) throw ()
  {
    return period_index_[axis] > 0;
  }

private:

  void init_grid_parents_() throw();
  void init_grid_levels_() throw();
  void init_grid_children_() throw();
  void init_grid_neighbors_() throw();
  void init_grid_faces_ (Domain & domain, Mpi & mpi) throw();
  void geomview_grid_ (FILE *fpr, bool full=true) throw ();

  /// Initialize il0_[], n0_[], and period_index[]
  void init_indices_ (bool is_periodic) throw();

  /// Initialize xl_[], xu_[], and period_domain[]
  void init_extents_ (bool is_periodic) throw();

  void insert_in_level_ (int level, Grid & grid) throw ();

  void deallocate_() throw();

};

//----------------------------------------------------------------------

class ItHierarchyGridsLocal
{

  /// @class ItHierarchyGridsLocal
  /// @date 2007-05-02
  /// @brief Iterator class for visiting all processor-local Grids in a Hierarchy

private:

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  ItHierarchyGridsLocal (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyGridsLocal (const ItHierarchyGridsLocal & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}
  void operator=(const ItHierarchyGridsLocal & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyGridsLocal () throw () {};
  
  /// Iterate through local Grids in the Hierarchy.
  Grid * operator++ (int) { 

    if (curr_ == hierarchy_->grids0_.size()) curr_ = 0;
    while (hierarchy_->grids0_[curr_] && ! hierarchy_->grids0_[curr_]->is_local()) curr_++;
    curr_ ++;
    return hierarchy_->grids0_[curr_-1];
  }

};

//----------------------------------------------------------------------

class ItHierarchyGridsAll
{

  /// @class ItHierarchyGridsAll
  /// @date 2007-05-02
  /// @brief Iterator class for visiting all Grids, both processor-local and remote, in a Hierarchy

private:

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  ItHierarchyGridsAll (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyGridsAll (const ItHierarchyGridsAll & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}
  void operator=(const ItHierarchyGridsAll & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyGridsAll () throw () {};
  
  /// Iterate through all Grids in the Hierarchy.
  Grid * operator++ (int) { 

    if (curr_ == hierarchy_->grids0_.size()) curr_ = 0;
    curr_ ++;
    return hierarchy_->grids0_[curr_-1];
  }

};

//----------------------------------------------------------------------


class ItHierarchyLevels
{

  /// @class ItHierarchyLevels
  /// @date 2007-05-02
  /// @brief Iterator class for visiting all Levels in a Hierarchy from coarsest to finest

private:

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  ItHierarchyLevels (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyLevels (const ItHierarchyLevels & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}

  void operator=(const ItHierarchyLevels & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyLevels () throw () {};
  
  /// Iterate through all Levels in the Hierarchy.
  Level * operator++ (int) { 

    if (curr_ == hierarchy_->levels0_.size()) curr_ = 0;
    curr_ ++;
    return hierarchy_->levels0_[curr_-1];
  }

};
#endif /* NEWGRAV_HIERRACHY_H */
