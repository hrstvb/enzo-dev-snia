// $Id: newgrav-domain.h 19 2010-03-19 23:04:28Z jbordner $

/// @file      newgrav-domain.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Declaration of the Domain class

#ifndef NEWGRAV_DOMAIN_H
#define NEWGRAV_DOMAIN_H

class Domain 
{

  /// @class    Domain
  /// @brief    Class for representing a problem domain's dimensionality and extents

private:

  int d_;
  Scalar xl_[3];
  Scalar xu_[3];

public:

  Domain () throw ();
  Domain (std::string parms) throw ();
  Domain (int d, Scalar xl[3], Scalar xu[3]) throw ();
  ~Domain () throw ();

  void print () throw ();
  void write (FILE *fp = 0) throw ();
  void input (std::string parms) throw ();

  void set_lower (Scalar xl_0, Scalar xl_1, Scalar xl_2) 
  { xl_[0] = xl_0;
    xl_[1] = xl_1;
    xl_[2] = xl_2; };

  void set_upper (Scalar xu_0, Scalar xu_1, Scalar xu_2) 
  { xu_[0] = xu_0;
    xu_[1] = xu_1;
    xu_[2] = xu_2; };

  /// Return the coordinates of the lower grid vertex.  No error checking on i.
  void lower(Scalar &xl_0, Scalar &xl_1, Scalar &xl_2) const throw ()
  { xl_0 = xl_[0];
    xl_1 = xl_[1];
    xl_2 = xl_[2]; };


  /// Return the coordinates of the upper grid vertex.  No error checking on i.
  void upper(Scalar &xu_0, Scalar &xu_1, Scalar &xu_2) const throw ()
  { xu_0 = xu_[0];
    xu_1 = xu_[1];
    xu_2 = xu_[2]; };



};

#endif /* NEWGRAV_DOMAIN_H */
