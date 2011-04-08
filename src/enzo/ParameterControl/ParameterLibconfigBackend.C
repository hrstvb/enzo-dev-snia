/* libconfig style parameter file wrapper */

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "ParameterBackends.h"

#include "libconfig/libconfig.h"
//using namespace libconfig;

class enzo_libconfig_backend : public interpreter
{
protected:
	std::string fname_;
	config_t cfg_;
	
	char strbuf_[2048];
	
	/*void get_path( const std::string& instr, std::string& path, std::string& key )
	{
		size_t i=instr.find_last_of('.');
		path = instr.substr(0,i);
		key  = instr.substr(i+1);
	}*/
	
	/*std::string conv( const Setting& elem )
	{
		int intval; double doubleval; std::string stringval;
		std::stringstream returnval;
		
		switch( elem.getType() )
		{
			case Setting::TypeInt: intval = elem; returnval << intval; break;
			case Setting::TypeFloat: doubleval = elem; returnval << std::setw(20) << std::setprecision(19) << doubleval; break;
			case Setting::TypeString: std::string dummy = elem; return dummy; break;
		}
		
		return returnval.str();
	}*/
	
	
public:
	explicit enzo_libconfig_backend( std::string fname )
	: interpreter(fname), fname_(fname)
	{
		config_init(&cfg_);
		
		if(! config_read_file(&cfg_, fname_.c_str()))
		{
			fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg_),
					config_error_line(&cfg_), config_error_text(&cfg_));
			throw std::runtime_error("parse error");
		}
		
		/*try{
			//cfg_.readFile( fname_.c_str() );
			//cfg_.setAutoConvert( true );	
			
		}
		catch( ParseException p)
		{
			std::cerr << "enzo_libconfig_backend: ";
			std::cerr << p.getError() << " " << p.getFile() << ":" << p.getLine() << std::endl;
		}*/
		
	}
	
	
	~enzo_libconfig_backend()
	{ 
		config_destroy(&cfg_);
	}
	
	void query( const std::string key, std::string &ret )
	{
		std::cerr << "accessing scalar " << key << std::endl;
		
		int intval; double doubleval; char* stringval;
		std::stringstream returnval;
		
		config_setting_t *setting = config_lookup(&cfg_, key.c_str() );
		
		if( setting == NULL )
		{	
			fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
			throw std::runtime_error("element not found");
		}
		
		switch( config_setting_type(setting) )
		{
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
				
		}
		//returnval << "1";
		
	}
#if 0
	void query( const std::string key, std::string& ret )
	{   
		/*const char *buf;
		if( !config_lookup_string(&cfg_, key.c_str(), &buf) )
		{	
			fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
			throw std::runtime_error("element not found");
		}
		ret = buf;*/
		std::cerr << key << std::endl;
		ret = lookup( key );
		
		/*
		try{
			const Setting &set = cfg_.getRoot()[key];
			ret = conv(set);			
			std::cerr << "successful access to \'" << key << "\'\n";
		}catch( SettingNotFoundException p )
		{
			std::cerr << "enzo_libconfig_backend: ";
			std::cerr << p.what() << " when looking up \'" << key << "\'" << std::endl;
		}*/
		

		//cfg_.lookupValue(key.c_str(),ret); 
	}
#endif
	
	void query_list( const std::string key, std::vector< std::string >& ret )
	{
		
		std::cerr << "accessing array " << key << std::endl;
		config_setting_t *setting;
		setting = config_lookup(&cfg_, key.c_str());
		
		if(setting != NULL)
		{
			int count = config_setting_length(setting);
			for( int i=0; i<count; ++i )
			{
				int intval; double doubleval; char* stringval;
				std::stringstream returnval;
				
				switch( config_setting_type( config_setting_get_elem(setting,i) ) )
				{
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
				}
			}
		}else{	
			fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
			throw std::runtime_error("element not found");
		}
		
			
		
		/*try{
			const Setting &array = cfg_.getRoot()[key];
			
			int count = array.getLength();
			
			for( int i=0; i<count; ++i )
			{
				const Setting &elem = array[i];
				std::string tmp = conv(elem);
				ret.push_back(tmp);
			}
		}catch( SettingNotFoundException p )
		{
			std::cerr << "enzo_libconfig_backend: ";
			std::cerr << p.what() << " when looking up \'" << key << "\'" << std::endl;
		}
		*/
		
		
	}
	
	size_t size( const std::string key )
	{
		/*Setting &array;
		size_t len = 0;
		try{
			std::string element_key;
			array = open_group( key, element_key );
			len = array[element_key].getLength();//array.getLength();
			std::cerr << "queried length of \'" << key << "\' is " << len << std::endl;
		}catch( SettingNotFoundException p )
		{
			std::cerr << "enzo_libconfig_backend: ";
			std::cerr << p.what() << " when looking up \'" << key << "\'" << std::endl;
		}
		return len;*/
		
		
		config_setting_t *setting;
		setting = config_lookup(&cfg_, key.c_str());
		int len = 0;
		if(setting != NULL)
		{
			len = config_setting_length(setting);
		}else{
			fprintf(stderr, "No '%s' setting in configuration file.\n",key.c_str());
			throw std::runtime_error("element not found");
		}
			
		std::cerr << "length of " << key << " = " << len;
		return len;
	}
};


namespace{
	interpreter_creator_concrete< enzo_libconfig_backend > c00("enzo2_libconfig");
}
