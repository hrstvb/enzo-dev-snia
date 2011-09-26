/* libconfig style parameter file wrapper */

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
  
  int update_parameters(config_setting_t *setting, char *path) {
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
    
    int type = setting->type;
    config_value_t value = setting->value;
    config_list_t *list = value.list;
      

  // ***** If this setting is a group, then recurse further down. *****

    if(type == CONFIG_TYPE_GROUP) {
      //    printf("%s is a group.\n",this_path);      
      
      if(list) {
	int len = list->length;	
	config_setting_t **s;
	for(s = list->elements; len--; s++)
	  update_parameters(*s, this_path);  // recursion
      }
      
      return 1;
    }
    

  // ***** If we get here then we're a leaf, i.e. a parameter. ********

    config_setting_t *old_setting;


    /* Check to see if this parameter already exists. (optional) */
    old_setting = config_lookup(&cfg_, this_path);

    if(old_setting) {
      printf("Parameter \"%s\" already exists -- updating value.\n",this_path);
    } else {
      printf("Adding new parameter \"%s\".\n",this_path);    
    }
    

    // Create a std::string version of the full path name.
    std::string setting_name = std::string( this_path );


    // Update the setting, or create it if necessary.
    switch( type ) {
    case CONFIG_TYPE_INT:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_INT );
      config_setting_set_int( old_setting, value.ival );
      break;
    case CONFIG_TYPE_FLOAT:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_FLOAT );
      config_setting_set_float( old_setting, value.fval );
      break;
    case CONFIG_TYPE_STRING:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_STRING );
      config_setting_set_string( old_setting, value.sval );
      break;
    case CONFIG_TYPE_BOOL:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_BOOL );
      config_setting_set_bool( old_setting, value.ival );
      break;
    case CONFIG_TYPE_ARRAY:
      old_setting = create_branch_if_necessary( setting_name, CONFIG_TYPE_ARRAY );
      
      // Remove all existing array elements from the old setting.
      while( config_setting_remove_elem(old_setting, 0) );
      
      // If there is more than one element...
      if(list) {
	config_setting_t *element;
	config_setting_t **s;
	int len = list->length;
	
	for(s = list->elements; len--; s++) {
	  
	  switch( (*s)->type ) {
	  case CONFIG_TYPE_INT:
	    element = config_setting_add(old_setting, NULL, CONFIG_TYPE_INT);
	    config_setting_set_int( element, ((*s)->value).ival );
	    break;
	  case CONFIG_TYPE_FLOAT:
	    element = config_setting_add(old_setting, NULL, CONFIG_TYPE_FLOAT );
	    config_setting_set_float( element, ((*s)->value).fval );
	    break;
	  case CONFIG_TYPE_STRING:
	    element = config_setting_add(old_setting, NULL, CONFIG_TYPE_STRING );
	    config_setting_set_string( element, ((*s)->value).sval );
	    break;
	  case CONFIG_TYPE_BOOL:
	    element = config_setting_add(old_setting, NULL, CONFIG_TYPE_BOOL );
	    config_setting_set_bool( element, ((*s)->value).ival );
	    break;
	  default:
	    fprintf(stderr,"unknown type: %d\n", (*s)->type);
	  }

	} // loop over all (new) array elements.
      }
      
      break;
    default:
      fprintf(stderr,"unknown type: %d\n", type);
    }
    
    
    return 1;
  }
  
  
public:
  
  explicit enzo_libconfig_backend( std::string fname, std::string defaults )
    : interpreter(fname), fname_(fname)
  {
    config_init(&cfg_);
    
    // We first "update" the newly created config object with the
    // defaults string. This is done first, so the user specified
    // parameters can overwrite these defaults.
    update(defaults);
    
    // Now we update with the user-specified parameters.
    update(fname_, true);    
  }
  
  
  ~enzo_libconfig_backend()
  { 
    config_destroy(&cfg_);
  }
  
  
  // This updates the existing config object with new values from
  // <input_string>, which must be in libconfig format.
  // 
  // When from_file==true, then <input_string> is treated as a
  // filename, from which to read the new values.
  //
  bool update( std::string input_string, bool from_file=false )
  {	
    
    // First we initialize a new (temporary) config
    // object.
    config_t tmp_cfg_;
    config_init(&tmp_cfg_);
    
    // Now we read in the new parameters into the temporary config.
    if( from_file ) {
      if(! config_read_file(&tmp_cfg_, input_string.c_str())) {
	fprintf(stderr, "%s:%d - %s\n", config_error_file(&tmp_cfg_), config_error_line(&tmp_cfg_), config_error_text(&tmp_cfg_));
	throw std::runtime_error("parse error");
      }

      if(! config_read_string(&tmp_cfg_, input_string.c_str())) {
	fprintf(stderr, "string:%d - %s\n", config_error_line(&tmp_cfg_), config_error_text(&tmp_cfg_));
	throw std::runtime_error("parse error");
      }
    }
    
    // Next we update the parameters in the original config object
    // with the new ones from the temporary one.
    update_parameters(tmp_cfg_.root, NULL);
        
    // Destroy temporary config object.
    config_destroy(&tmp_cfg_);
    
    return true;
  }	
  
  bool dump( std::string fname )  
  {
    if(! config_write_file( &cfg_, fname.c_str() ) )
      return false;
    
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
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_INT );
    return config_setting_set_int( setting, value )==CONFIG_TRUE;
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
  
  bool set( std::string key, std::string value )
  {
    config_setting_t *setting = create_branch_if_necessary( key, CONFIG_TYPE_STRING );
    return config_setting_set_string( setting, value.c_str() )==CONFIG_TRUE;
  }
  
  bool set_list( std::string key, size_t n, const int* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *setting = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      setting = config_setting_add(array, NULL, CONFIG_TYPE_INT);
      ret &= config_setting_set_int(setting, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const long long* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *setting = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      setting = config_setting_add(array, NULL, CONFIG_TYPE_INT64);
      ret &= config_setting_set_int64(setting, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const bool* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *setting = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      setting = config_setting_add(array, NULL, CONFIG_TYPE_BOOL);
      ret &= config_setting_set_int64(setting, (int)values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const double* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *setting = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      setting = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
      ret &= config_setting_set_float(setting, values[i])==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  bool set_list( std::string key, size_t n, const std::string* values )
  {
    bool ret = true;
    
    config_setting_t *array = create_branch_if_necessary( key, CONFIG_TYPE_ARRAY );
    config_setting_t *setting = NULL;
    
    for( size_t i=0l; i<n; ++i ) {
      setting = config_setting_add(array, NULL, CONFIG_TYPE_STRING);
      ret &= config_setting_set_string(setting, values[i].c_str())==CONFIG_TRUE;
    }
    
    return ret;
  }
  
  int query( std::string key, std::string &ret )
  {
    std::cerr << "accessing scalar " << key << std::endl;
    
    int intval; double doubleval; char* stringval;
    std::stringstream returnval;
    
    config_setting_t *setting = config_lookup(&cfg_, key.c_str() );
    
    if( setting == NULL ) {	
      //		  fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
      //		  throw std::runtime_error("element not found");
      return 0;
    }
    
    switch( config_setting_type(setting) ) {
    case CONFIG_TYPE_INT: 
      intval = config_setting_get_int(setting); 
      returnval << intval; 
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
    
    std::cerr << "accessing array " << key << std::endl;
    config_setting_t *setting;
    setting = config_lookup(&cfg_, key.c_str());
    
    if(setting != NULL) {
      int count = config_setting_length(setting);
      for( int i=0; i<count; ++i ) {
	int intval; double doubleval; char* stringval;
	std::stringstream returnval;
	
	switch( config_setting_type( config_setting_get_elem(setting,i) ) ) {
	case CONFIG_TYPE_INT: 
	  intval = config_setting_get_int_elem(setting,i); 
	  returnval << intval; 
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
      //		  fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
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
    
    std::cerr << "length of " << key << " = " << len;
    return len;
  }
};


namespace {
  interpreter_creator_concrete< enzo_libconfig_backend > c00("enzo2_libconfig");
}
