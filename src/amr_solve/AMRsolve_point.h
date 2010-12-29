/// @file      newgrav-point.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Point class

#ifndef AMRSOLVE_POINT_H
#define AMRSOLVE_POINT_H

class AMRsolve_Grid;

class AMRsolve_Point 
{

  /// @class    AMRsolve_Point
  /// @brief    Declaration of AMRsolve_Point class for representing a point mass

private:

  // AMRsolve_Point attributes known when created

  static int d_;

  Scalar m_;
  Scalar *x_;
  int igrid_;

public:

  AMRsolve_Point(Scalar m=0.0, Scalar x1=0.0, Scalar x2=0.0, Scalar x3=0.0) throw();
  AMRsolve_Point(std::string parms) throw();
  ~AMRsolve_Point() throw();

  void print() throw();
  void write(FILE* fp=0) throw();
  void input(std::string parms) throw();

  /// Set the dimension of the space that the AMRsolve_Point is contained in
  static void set_dim(int d) { d_ = d; }

  /// Return the id of the containing Grid
  int igrid() {return igrid_; };

  /// Return the mass of the point
  Scalar mass() {return m_; };

  /// Return the ith coordinate of the point.  No error checking.
  Scalar x(int i) {return x_[i]; };

private:

  void alloc_() throw();
  void dealloc_() throw();

};

#endif /* AMRSOLVE_POINT_H */
