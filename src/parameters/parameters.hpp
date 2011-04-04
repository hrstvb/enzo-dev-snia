// $Id: parameters.hpp 2044 2011-03-02 19:35:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

/// @file     parameters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:44:21 PDT 2009
/// @brief    Include file for the \ref Parameters component

/// @enum     parameter_enum
/// @brief    Parameter data type
enum parameter_enum {
  parameter_unknown,
  parameter_integer,
  parameter_scalar,
  parameter_string,
  parameter_logical,
  parameter_list,
};

extern const char * parameter_type_name [];

class ExceptionParametersBadType {};


//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "parse.h"
#include "parameters_Param.hpp"
#include "parameters_ParamNode.hpp"
#include "parameters_Parameters.hpp"

#endif /* PARAMETERS_HPP */

