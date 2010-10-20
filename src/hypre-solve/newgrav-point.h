// $Id: newgrav-point.h 19 2010-03-19 23:04:28Z jbordner $

/// @file      newgrav-point.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Declaration of the Point class

#ifndef NEWGRAV_POINT_H
#define NEWGRAV_POINT_H

class Grid;

class Point 
{

  /// @class    Point
  /// @brief    Declaration of Point class for representing a point mass

private:

  // Point attributes known when created

  static int d_;

  Scalar m_;
  Scalar *x_;
  int igrid_;

public:

  Point (Scalar m=0.0, Scalar x1=0.0, Scalar x2=0.0, Scalar x3=0.0) throw () ;
  Point (std::string parms) throw ();
  ~Point () throw ();

  void print () throw ();
  void write (FILE *fp = 0) throw ();
  void input (std::string parms) throw ();

  /// Set the dimension of the space that the Point is contained in
  static void set_dim (int d) 
  { d_ = d; }

  /// Return the id of the containing Grid
  int igrid () 
  {return igrid_; } ;

  /// Return the mass of the point
  Scalar mass () 
  {return m_; } ;

  /// Return the ith coordinate of the point.  No error checking.
  Scalar x (int i) 
  {return x_[i]; } ;

private:

  void alloc_ () throw ();
  void dealloc_ () throw ();

};

#endif /* NEWGRAV_POINT_H */
