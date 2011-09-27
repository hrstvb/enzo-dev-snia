// libconfig style parameter file wrapper

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "string.h"

#include "ParameterBackends.h"

#include "libconfig/libconfig.h"
//using namespace libconfig;

// set this for verbose comments to stdout.
int verbose_paramconfig = 0;

const char* var_types[9] = { "NONE", "GROUP", "INT", "INT64", "FLOAT", "STRING", "BOOL", "ARRAY", "LIST" };

class enzo_libconfig_backend : public interpreter
{
protected:
  std::string fname_;
  config_t cfg_;
  
  char strbuf_[2048];
  
  
protected:
  
  // splits input string 'str' into substrings 'results' delimited by 'delim'
  size_t string_split( std::string str, const std::string delim, std::vector<std::string>& results ) {
    size_t i0 = 0;
    while( (i0 = str.find_first_of(delim)) != str.npos ) {
      if( i0 > 0l )
	results.push_back( str.substr(0,i0) );
      str = str.substr(i0+1);
    }
    
    if( str.length() > 0l )
      results.push_back(str);
    
    return results.size();
  }
  
  // Checks that the types between two setting match.
  int check_types(config_setting_t *s_old, config_setting_t *s_new, char *path) {
    short old_type = s_old->type;
    short new_type = s_new->type;
    
    // A type mismatch throws an error, unless it's one of these two
    // cases:
    //
    // (a) old_type == float, new_type == int64
    // (b) old_type == bool, new_type == int64
    
    if(old_type != new_type) {

      if( (old_type == CONFIG_TYPE_FLOAT) && 
	  (new_type == CONFIG_TYPE_INT64) ) {

	s_new->type = old_type;
	(s_new->value).fval = (double) (s_new->value).llval;
      } else if( (old_type == CONFIG_TYPE_BOOL) &&
		 (new_type == CONFIG_TYPE_INT64) ) {

	s_new->type = old_type;
      } else {
	std::string error_message;
	
	error_message = "Type mismatch for parameter \"" + std::string(path) + "\": old type = " + std::string(var_types[old_type]) + ", new type = " + std::string(var_types[new_type]) + ".";
	
	throw std::runtime_error(error_message);
      }
    }

    return 1;
  }


  // returns the setting pointer to the element, creating all groups along the way
  config_setting_t *create_branch_if_necessary( std::string key, int setting_type ) {
    std::vector<std::string> split_key;
    config_setting_t *sptr = config_root_setting( &cfg_ );
    config_setting_t *next = NULL;
    
    string_split( key, ".", split_key );
    
    for( size_t i=0; i<split_key.size()-1; ++i ) {
      next = config_setting_get_member( sptr, split_key[i].c_str() );
      if( next == NULL )
	next = config_setting_add( sptr, split_key[i].c_str(), CONFIG_TYPE_GROUP );
      sptr = next;
    }
    
    next = config_setting_get_member( sptr, split_key.back().c_str() );
    if( next == NULL )
      next = config_setting_add( sptr, split_key.back().c_str(), setting_type );
    
    return next;
  }
  
  int update_parameters(config_setting_t *setting, char *path, bool overwrite=false) {
    char this_path[MAX_PARAM_LENGTH+1];
    this_path[0] = 0;
    
    if(setting->name) {
      if(path[0]==0) {
	
	//make sure new path length does not exceed MAX_PARAM_LENGTH
	if(strlen(setting->name) > MAX_PARAM_LENGTH)
	  throw std::runtime_error("Path length exceeds MAX_PARAM_LENGTH!\n");
	
	sprintf(this_path,"%s",setting->name);
      }
      else {
	
	//make sure new path length does not exceed MAX_PARAM_LENGTH
	if(strlen(path) + strlen(setting->name) > MAX_PARAM_LENGTH)
	  throw std::runtime_error("Path length exceeds MAX_PARAM_LENGTH!\n");
	
	sprintf(this_path,"%s.%s",path,setting->name);
      }
    }
    
    short *type = &(setting->type);
    config_value_t *value = &(setting->value);
    config_list_t *list = value->list;
      

  // ***** If this setting is a group, then recurse further down. *****

    if((*type) == CONFIG_TYPE_GROUP) {
      if(verbose_paramconfig)
	printf("update_parameters: %s is a group.\n",this_path);      
      
      if(list) {
	int len = list->length;	
	config_setting_t **s;
	for(s = list->elements; len--; s++)
	  update_parameters(*s, this_path, overwrite);  // recursion
      }
      
      return 1;
    }
    

  // ***** If we get here then we're a leaf, i.e. a parameter. ********

    // Force integers to be 64bit.
    if( (*type) == CONFIG_TYPE_INT ) {
      (*type) = CONFIG_TYPE_INT64;
      value->llval = (long long) (value->ival);
    }

    // Check to see if this parameter already exists.
    config_setting_t *old_setting = config_lookup(&cfg_, this_path);

    if( (old_setting) && (! overwrite) ) {
      if(verbose_paramconfig)
	printf("Parameter \"%s\" already exists -- skipping it.\n",this_path);
      return 1;
    }

    if(old_setting) {
      if(verbose_paramconfig)
	printf("Parameter \"%s\" already exists -- updating value.\n",this_path);
      check_types(old_setting, setting, this_path);
    } else {
      if(verbose_paramconfig)
	printf("Adding new parameter \"%s\".\n",this_path);    
    }
        
    // Create a std::string version of the full path name.
    std::string setting_name = std::string( this_path );


    // Update the setting, or create it if necessary.
    switch( (*type) ) {
    case CONFIG_TYPE_INT:
      // should never get here, since we force 64bit ints above...
      throw std::runtime_error("32bit int detected! This should not have happened...");
      break;
    case CONFIG_TYPE_INT64:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_INT64 );
      config_setting_set_int64( old_setting, value->llval );
      break;
    case CONFIG_TYPE_FLOAT:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_FLOAT );
      config_setting_set_float( old_setting, value->fval );
      break;
    case CONFIG_TYPE_STRING:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_STRING );
      config_setting_set_string( old_setting, value->sval );
      break;
    case CONFIG_TYPE_BOOL:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_BOOL );
      config_setting_set_bool( old_setting, (int) value->llval );
      break;
    case CONFIG_TYPE_ARRAY:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_ARRAY );
      
      // Remove all existing array elements from the old setting.
      //
      // No. We cannot do that, since we need to preserve the full set
      // of default values for arrays like MinimumOverDensityForRefinement.

      // while( config_setting_remove_elem(old_setting, 0) );

      
      // If there is more than zero elements...
      if(list) {
	int Nelem = list->length;		

	config_setting_t **s;
	s = list->elements;

	short *this_type;	
	config_value_t *this_value;
	for(int i = 0; i < Nelem; i++, s++) {

	  this_value = &((*s)->value);

	  // Force integers to be 64bit. (This needs to happen inside
	  // this loop in order to typecast the value to 'long long'.)
	  this_type = &((*s)->type);
	  if( (*this_type) == CONFIG_TYPE_INT) {
	    (*this_type) = CONFIG_TYPE_INT64;
	    this_value->llval = (long long) this_value->ival;
	  }

	  config_setting_t *element;
	  element = config_setting_get_elem(old_setting, i);

	  if(element)
	    check_types(element, (*s), this_path);

	  switch( (*this_type) ) {
	  case CONFIG_TYPE_INT:
	    // should never get here, since we force 64bit ints above...
	    throw std::runtime_error("32bit int detected! This should not have happened...");
	    break;
	  case CONFIG_TYPE_INT64:
	    if(!element)
	      element = config_setting_add(old_setting, NULL, CONFIG_TYPE_INT64);
	    config_setting_set_int64( element, this_value->llval );
	    break;
	  case CONFIG_TYPE_FLOAT:
	    if(!element)
	      element = config_setting_add(old_setting, NULL, CONFIG_TYPE_FLOAT );
	    config_setting_set_float( element, this_value->fval );
	    break;
	  case CONFIG_TYPE_STRING:
	    if(!element)
	      element = config_setting_add(old_setting, NULL, CONFIG_TYPE_STRING );
	    config_setting_set_string( element, this_value->sval );
	    break;
	  case CONFIG_TYPE_BOOL:
	    if(!element)
	      element = config_setting_add(old_setting, NULL, CONFIG_TYPE_BOOL );
	    config_setting_set_bool( element, (int) this_value->llval );
	    break;
	  default:
	    std::string error_message;
	    error_message = "Unknown type for parameter \"" + std::string(this_path) + "\": " + std::string(var_types[(*this_type)]) + ".";
	    throw std::runtime_error(error_message);
	  }

	} // loop over all (new) array elements.
      }
      
      break;
    default:
      std::string error_message;
      error_message = "Unknown type for parameter \"" + std::string(this_path) + "\": " + std::string(var_types[(*type)]) + ".";
      throw std::runtime_error(error_message);
    }
    
    
    return 1;
  }
  
  
public:
  
  explicit enzo_libconfig_backend( std::string fname, std::string defaults )
    : interpreter(fname), fname_(fname) {
    if(verbose_paramconfig)
      printf("Initializing parameter config (libconfig backend).\n");

    config_init(&cfg_);
    
    // We first "update" the newly created config object with the
    // defaults string. This is done first in order to set some arrays
    // (like CellFlaggingMethod, etc.) which require more elements
    // than the user typically specifies.
    update(defaults);
    
    // useful for debugging:
    //    dump("dump_before.cfg");

    // Now we update (from_file=true, overwrite=true) with the
    // user-specified parameters.
    update(fname_, true, true);    

    // useful for debugging:
    //    dump("dump.cfg");
  }
  
  
  ~enzo_libconfig_backend()
  { 
    config_destroy(&cfg_);
  }
  
  
  // This updates the existing config object with new values from
  // <input_string>, which must be in libconfig format.
  //
  // Unless overwrite=true, this will simply skip parameters that have
  // already been defined. This is so UpdateDefaults() doesn't wipe
  // out the user's parameter choices.
  // 
  // When from_file==true, then <input_string> is treated as a
  // filename, from which to read the new values.
  //
  bool update( std::string input_string, bool from_file=false, bool overwrite=false )
  {	
    
    // First we initialize a new (temporary) config
    // object.
    config_t tmp_cfg_;
    config_init(&tmp_cfg_);
    
    // Now we read in the new parameters into the temporary config.
    if( from_file ) {

      if(verbose_paramconfig)
	printf("Updating parameters from file <%s>.\n", input_string.c_str());

      if(! config_read_file(&tmp_cfg_, input_string.c_str())) {
	fprintf(stderr, "%s:%d - %s\n", config_error_file(&tmp_cfg_), config_error_line(&tmp_cfg_), config_error_text(&tmp_cfg_));
	throw std::runtime_error("parse error");
      }
    } else {

      if(verbose_paramconfig)
	printf("Updating parameters from string:\n\"%s\n\".\n", input_string.c_str());

      if(! config_read_string(&tmp_cfg_, input_string.c_str())) {
	fprintf(stderr, "string:%d - %s\n", config_error_line(&tmp_cfg_), config_error_text(&tmp_cfg_));
	throw std::runtime_error("parse error");
      }
    }
    
    // Next we update the parameters in the original config object
    // with the new ones from the temporary one.
    update_parameters(tmp_cfg_.root, NULL, overwrite);
        
    // Destroy temporary config object.
    config_destroy(&tmp_cfg_);
    
    return true;
  }	
  
  bool dump( char fname[], char header_string[]=NULL )
  {
    if(! config_write_file( &cfg_, fname ) )
      return false;
    
    // prepend a header string, if specified.
    if(header_string) {
      FILE *fp;

      // find number of characters
      char c;
      int Nchar = 0;
      fp = fopen(fname, "rt");      
      while((c = fgetc(fp)) != EOF) Nchar++;
      fclose(fp);

      // read the entire file into a buffer
      char *buffer = new char[Nchar];
      fp = fopen(fname, "rt");      
      for(int i=0;i<Nchar;i++)
	buffer[i] = fgetc(fp);
      fclose(fp);
      
      // write out the header and then the buffer
      fp = fopen(fname, "wt");      
      fprintf(fp,"%s\n",header_string);
      for(int i=0;i<Nchar;i++)
	fputc(buffer[i],fp);
      fclose(fp);
      
      delete [] buffer;
    }

    return true;
  }
  
  bool remove( std::string key )
  {
    config_setting_t *setting = config_lookup(&cfg_, key.c_str() );
    if( setting == NULL )
      return true;
    
    if( config_setting_remove( config_root_setting(&cfg_), key.c_str() ) == CONFIG_TRUE )
      return true;
    
    return false;
  }
  
  bool set( std::string key, int value )
  {

    // We force 64bit integers!
    long long lvalue = (long long) value;
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_INT64 );
    return config_setting_set_int64( setting, lvalue )==CONFIG_TRUE;
    
    // config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_INT );
    // return config_setting_set_int( setting, value )==CONFIG_TRUE;
  }
  
  bool set( std::string key, long long value )
  {
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_INT64 );
    return config_setting_set_int64( setting, value )==CONFIG_TRUE;
  }
  
  bool set( std::string key, bool value )
  {
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_BOOL );
    return config_setting_set_bool( setting, (int)value )==CONFIG_TRUE;
  }
  
  bool set( std::string key, double value )
  {
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_FLOAT );
    return config_setting_set_float( setting, value )==CONFIG_TRUE;
  }
  
  bool set( std::string key, char value[] )
  {

    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_STRING );
    return config_setting_set_string( setting, value )==CONFIG_TRUE;
  }
  
  bool set_list( std::string key, size_t n, const int* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *element = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      element = config_setting_get_elem(array, i);

      // We force 64bit integers!
      long long lvalue = (long long) values[i];
      if(!element)
	element = config_setting_add(array, NULL, CONFIG_TYPE_INT64);
      ret &= config_setting_set_int64(element, lvalue)==CONFIG_TRUE;

      // element = config_setting_add(array, NULL, CONFIG_TYPE_INT);
      // ret &= config_setting_set_int(element, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const long long* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *element = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      element = config_setting_get_elem(array, i);

      if(!element)
	element = config_setting_add(array, NULL, CONFIG_TYPE_INT64);
      ret &= config_setting_set_int64(element, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const bool* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *element = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      element = config_setting_get_elem(array, i);

      if(!element)
	element = config_setting_add(array, NULL, CONFIG_TYPE_BOOL);
      ret &= config_setting_set_int64(element, (int)values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const double* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *element = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      element = config_setting_get_elem(array, i);

      if(!element)
	element = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
      ret &= config_setting_set_float(element, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, char* values[] )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *element = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      element = config_setting_get_elem(array, i);

      if(!element)
	element = config_setting_add(array, NULL, CONFIG_TYPE_STRING);
      ret &= config_setting_set_string(element, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  int query( std::string key, std::string &ret )
  {
    if(verbose_paramconfig)
      printf("accessing scalar parameter \"%s\".\n", key.c_str());
    
    int intval; double doubleval; char* stringval;
    long long longval;
    std::stringstream returnval;
    
    config_setting_t *setting = config_lookup(&cfg_, key.c_str() );
    
    if( setting == NULL ) {	
      fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
      //		  throw std::runtime_error("element not found");
      return 0;
    }
    
    switch( config_setting_type(setting) ) {
    case CONFIG_TYPE_INT: 
      // should never get here, since we force 64bit ints...
      throw std::runtime_error("32bit int detected! This should not have happened...");

      // intval = config_setting_get_int(setting); 
      // returnval << intval; 
      // ret =  returnval.str();
      break;
    case CONFIG_TYPE_INT64: 
      longval = config_setting_get_int64(setting); 
      returnval << longval; 
      ret =  returnval.str();
      break;
    case CONFIG_TYPE_FLOAT: 
      doubleval = config_setting_get_float(setting); 
      returnval << std::setw(20) << std::setprecision(19) << doubleval; 
      ret =  returnval.str();
      break;
    case CONFIG_TYPE_STRING: 
      ret = config_setting_get_string(setting); 
      break;
    case CONFIG_TYPE_BOOL:
      if( config_setting_get_bool(setting) )
	ret = "1";
      else
	ret = "0";
      break;
    default:
      std::cerr << "unknown type" << config_setting_get_format(setting) << "\n";
      return 0;
      
    }
    return 1;		
  }
  
  int query_list( std::string key, std::vector< std::string >& ret )
  {
    if(verbose_paramconfig)
      printf("accessing array parameter \"%s\".\n", key.c_str());
    config_setting_t *setting;
    setting = config_lookup(&cfg_, key.c_str());
    
    if(setting != NULL) {
      int count = config_setting_length(setting);
      for( int i=0; i<count; ++i ) {
	int intval; double doubleval; char* stringval;
	long long longval;
	std::stringstream returnval;
	
	switch( config_setting_type( config_setting_get_elem(setting,i) ) ) {
	case CONFIG_TYPE_INT: 
	  // should never get here, since we force 64bit ints...
	  throw std::runtime_error("32bit int detected! This should not have happened...");

	  // intval = config_setting_get_int_elem(setting,i); 
	  // returnval << intval; 
	  // ret.push_back(returnval.str());
	  break;
	case CONFIG_TYPE_INT64: 
	  longval = config_setting_get_int64_elem(setting,i); 
	  returnval << longval; 
	  ret.push_back(returnval.str());
	  break;
	case CONFIG_TYPE_FLOAT: 
	  doubleval = config_setting_get_float_elem(setting,i); 
	  returnval << std::setw(20) << std::setprecision(19) << doubleval;
	  ret.push_back(returnval.str());
	  break;
	case CONFIG_TYPE_STRING: 
	  ret.push_back( config_setting_get_string_elem(setting,i) ); 
	  break;
	case CONFIG_TYPE_BOOL:
	  if( config_setting_get_bool_elem(setting,i) )
	    ret.push_back("1");
	  else
	    ret.push_back("0");
	  break;
	default:
	  std::cerr << "unknown type\n";
	  return 0;
	}
      }
    } else {	
      fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
      //		  throw std::runtime_error("element not found");
      return 0;
    }
    return 1;
  }
  
  size_t size( std::string key ) {
    
    config_setting_t *setting;
    setting = config_lookup(&cfg_, key.c_str());
    int len = 0;
    if(setting != NULL) {
      len = config_setting_length(setting);
    } else {
      fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
      //		  throw std::runtime_error("element not found");
      return 0;
    }
        
    
    if(verbose_paramconfig)
      printf("Parameter \"%s\" has %d elements.\n", key.c_str(), len);
    return len;
  }
};


namespace {
  interpreter_creator_concrete< enzo_libconfig_backend > c00("enzo2_libconfig");
}
