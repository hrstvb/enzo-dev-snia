/// @file      AMRsolve_problem.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Problem class

#ifndef AMRGRAV_PROBLEM_H
#define AMRGRAV_PROBLEM_H

#define BUFFER_LENGTH 255

class AMRsolve_Problem {

  /// @class    AMRsolve_Problem
  /// @brief    Declaration of AMRsolve_Problem class for storing test problem parameters

  //----------------------------------------------------------------------


private:

  AMRsolve_Parameters           parameters_; // Raw input parameters
  std::vector<AMRsolve_Point *> points_;     // List of point masses
  AMRsolve_Hierarchy            hierarchy_;  // AMR mesh hierarchy
  AMRsolve_Domain               domain_;     // Problem domain

  //----------------------------------------------------------------------

public:

  AMRsolve_Problem() throw();

  ~AMRsolve_Problem() throw();

  AMRsolve_Problem(const AMRsolve_Problem&) throw();

  AMRsolve_Problem& operator=(const AMRsolve_Problem&) throw();

  //--------------------------------------------------

  void print() throw();

  void write(FILE *fp=0) throw();

  /// Read in a problem file
  void read(std::string filename) throw();

  /// Return the dimension
  int dimension() { return hierarchy_.dimension(); }

  /// Return the hierarchy
  AMRsolve_Hierarchy& hierarchy() { return hierarchy_; }

  /// Return the domain
  AMRsolve_Domain& domain() { return domain_; }

  /// Return a pointer to the ith AMRsolve_Point
  AMRsolve_Point& point(int i) { return *points_[i]; }

  /// Return vector of points
  std::vector<AMRsolve_Point *> points() { return points_; }

  /// Return the number of points
  int num_points() { return points_.size(); }

  /// Return a pointer to the ith Grid
  AMRsolve_Grid& return_grid(int i) { return hierarchy_.return_grid(i); }

  /// Return the number of Grids
  int num_grids() { return hierarchy_.num_grids(); }

  /// Return a pointer to the ith AMRsolve_Level
  AMRsolve_Level& level(int i) { return hierarchy_.level(i); }

  /// Return the number of Levels
  int num_levels() { return hierarchy_.num_levels(); }

  /// Return the AMRsolve_Parameters object
  AMRsolve_Parameters& parameters() { return parameters_; }

  void deallocate_() throw();
};

#endif /* AMRGRAV_PROBLEM_H */
