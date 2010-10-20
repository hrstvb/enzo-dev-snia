// $Id: newgrav-parameters.h 19 2010-03-19 23:04:28Z jbordner $

/// @file      newgrav-parameters.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Declaration of the Parameters class

#ifndef NEWGRAV_PARAMETERS_H
#define NEWGRAV_PARAMETERS_H

#define BUFFER_LENGTH 255

class Parameters 
{

  /// @class    Parameters
  /// @brief    Declaration of Parameters class for storing and accessing run-time parameters

  friend class ItParameters;

private:

  std::multimap<std::string, std::string> values_;

public:

  Parameters () throw () ;
  ~Parameters () throw () ;

  /// Read a file of key-value pairs.  May be more than one value per key.
  void read (std::string file) throw () ;

  /// Print all parameters to stdout.
  void print () throw () ;

  /// Associate the given value with the given key.
  void add_parameter (std::string key, std::string val) throw () ;

  /// Retrieve the ith value of the the given parameter.
  std::string ith_value  (std::string key, int i) const throw () ;

  /// Return the multiplicity of values for the given key.  May be 0.
  int num_values  (std::string key) const throw () ;

  /// Return the value for the given key.
  std::string value (std::string key) const throw ();

private:

  int readline_ (FILE *, char * buffer, int n) throw ();

};

//----------------------------------------------------------------------

class ItParameters
{
  /// @class ItParameters
  /// @date 2007-08-30
  /// @brief Iterate through each key-value pairs in a Parameters object
  ///
  /// Example usage:
  ///
  /// @verbatim
  ///
  /// ItParameters it_parameter (parameters);
  /// while (it_parameter++) {
  ///    std::string key   = it_parameter.key();
  ///    std::string value = it_parameter.value();
  ///    ... 
  /// }
  /// @endverbatim

private:

  Parameters * pparameters_;
  std::multimap<std::string, std::string>::iterator curr_;
  std::multimap<std::string, std::string>::iterator next_;

public:

  ItParameters (Parameters & parameters) throw ()
    : pparameters_(& parameters),
      curr_(),
      next_(parameters.values_.begin())
  { }

  ~ItParameters () throw () {};
  
  /// Iterate through all Grids in the Grid.
  int operator++ (int) { 
    
    curr_ = next_;
    int not_done = 1;
    if (next_ == pparameters_->values_.end()) {
      next_ = pparameters_->values_.begin();
      not_done = 0;
    }
    ++next_;
    return not_done;
  }

  std::string key () throw ()
  { return curr_->first; }
  std::string value () throw ()
  { return curr_->second; }

};

extern Parameters * hypre_parameters;

#endif /* NEWGRAV_PARAMETERS_H */
