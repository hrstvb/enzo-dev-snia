// $Id: parameters_Param.cpp 2169 2011-04-03 17:01:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Param.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sun Oct 11 15:02:08 PDT 2009
/// @brief    Implementation of the Param class

#include "parameters.hpp"

void Param::set (struct param_struct * node)
/// @param   node  The node from which to copy the type and value
{

  value_accessed_ = false;

  switch (node->type) {
  case enum_parameter_integer:
    set_integer_(node->integer_value);
    break;
  case enum_parameter_scalar:
    set_scalar_(node->scalar_value);
    break;
  case enum_parameter_string:
    set_string_(node->string_value);
    break;
  case enum_parameter_logical:
    set_logical_(node->logical_value);
    break;
  case enum_parameter_list:
    set_list_(node->list_value);
    break;
  case enum_parameter_unknown:
  case enum_parameter_sentinel:
  case enum_parameter_function:
  case enum_parameter_group:
  case enum_parameter_identifier:
    break;
  }
}

void Param::dealloc_() 
///
{ 
  switch (type_) {
  case parameter_string: 
    dealloc_string_(); 
    break;
  case parameter_list:   
    dealloc_list_(value_list_); 
    break;
  case parameter_unknown:
  case parameter_integer:
  case parameter_scalar:
  case parameter_logical:
    break;
  }
} 


void Param::write(FILE * file_pointer,
		  std::string parameter)
/// @param file_pointer  File pointer to which the parameter is written
/// @param parameter     Name of this parameter
/// @todo  Writing lists is not implemented yet
{

  // Write access indicator
  fprintf (file_pointer, value_accessed_ ? "[*] " : "[ ] ");

  // Write the parameter name
  fprintf (file_pointer,"%s ",parameter.c_str());

  // Write the parameter value
  fprintf (file_pointer,"%s", value_to_string().c_str());

  fprintf (file_pointer,"\n");
}


std::string Param::value_to_string ()
{
  char string_buffer[80];
  switch (type_) {
  case parameter_string: 
    sprintf (string_buffer,"%s",value_string_);
    break;
  case parameter_list:
    sprintf (string_buffer,"LIST\n");
    printf ("INCOMPLETE: Param::write");
    break;
  case parameter_integer:
    sprintf (string_buffer,"%d",value_integer_);
    break;
  case parameter_scalar:
    sprintf (string_buffer,"%g",value_scalar_);
    break;
  case parameter_logical:
    sprintf (string_buffer,"%s",value_logical_ ? "true" : "false");
    break;
  case parameter_unknown:
    sprintf (string_buffer,"UNKNOWN\n");
    break;
  }  
  return string_buffer;
}
 
void Param::dealloc_list_ (list_type * value)
/// @param value List to be deallocated
{
  for (unsigned i=0; i<(*value).size(); i++) {
    if ((*value)[i]->type_ == parameter_list) {
      dealloc_list_ ((*value)[i]->value_list_);
    } else {
      delete ( (*value)[i] );
    }
  }
  delete value;
}

