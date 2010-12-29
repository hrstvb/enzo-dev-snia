/// @file      newgrav-parameters.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Parameters class

#ifndef AMRSOLVE_PARAMETERS_H
#define AMRSOLVE_PARAMETERS_H

#define BUFFER_LENGTH 255

class AMRsolve_Parameters 
{

  /// @class    AMRsolve_Parameters
  /// @brief    Declaration of AMRsolve_Parameters class for storing and accessing run-time parameters

  friend class ItParameters;

private:

  std::multimap<std::string, std::string> values_;

public:

  AMRsolve_Parameters() throw();
  ~AMRsolve_Parameters() throw();

  /// Read a file of key-value pairs.  May be more than one value per key.
  void read(std::string file) throw();

  /// Print all parameters to stdout.
  void print() throw();

  /// Associate the given value with the given key.
  void add_parameter(std::string key, std::string val) throw();

  /// Retrieve the ith value of the the given parameter.
  std::string ith_value(std::string key, int i) const throw();

  /// Return the multiplicity of values for the given key.  May be 0.
  int num_values(std::string key) const throw();

  /// Return the value for the given key.
  std::string value(std::string key) const throw();

private:

  int readline_(FILE*, char* buffer, int n) throw();

};

//----------------------------------------------------------------------

class ItParameters
{
  /// @class ItParameters
  /// @date 2007-08-30
  /// @brief Iterate through each key-value pairs in a AMRsolve_Parameters object
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

  AMRsolve_Parameters* pparameters_;
  std::multimap<std::string, std::string>::iterator curr_;
  std::multimap<std::string, std::string>::iterator next_;

public:

  ItParameters(AMRsolve_Parameters& parameters) throw()
    : pparameters_(& parameters), curr_(), next_(parameters.values_.begin())
  { }

  ~ItParameters() throw() {};
  
  /// Iterate through all parameters in the Parameters.
  int operator++(int) { 
    curr_ = next_;
    int not_done = 1;
    if (next_ == pparameters_->values_.end()) {
      next_ = pparameters_->values_.begin();
      not_done = 0;
    }
    ++next_;
    return not_done;
  }

  std::string key() throw()
  { return curr_->first; }

  std::string value() throw()
  { return curr_->second; }

};

extern AMRsolve_Parameters* hypre_parameters;

#endif /* AMRSOLVE_PARAMETERS_H */
