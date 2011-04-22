#ifndef PARAMETER_CONTROL_H
#define PARAMETER_CONTROL_H

#include <map>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <cstdarg>
#include <cstdio>

#include <iostream>

#include "ParameterBackends.h"

/*********************************************************************/
// SOME HELPER FUNCTIONS

//! converts between different variable types
/*!
 *  The main purpose of this function is to parse and convert
 *  a string argument into numbers, booleans, etc...
 * @param ival the input value (typically a std::string)
 * @param oval the interpreted/converted value
 */
template <class in_value, class out_value>
inline void convert( const in_value & ival, out_value & oval) 
{
	std::stringstream ss;
	ss << ival; //.. insert value into stream
	ss >> oval; //.. retrieve value from stream
	
	if (! ss.eof()) {
		//.. conversion error
		//std::cerr << "Error: conversion of \'" << ival << "\' failed." << std::endl;
		throw std::runtime_error(std::string("invalid conversion to ")+typeid(out_value).name()+'.');
	}
}


template<>
inline void convert<std::string,std::string>( const std::string & ival, std::string & oval)
{
	oval = ival;
}



//... actual configuration class

class Configuration
{
	interpreter_creator *ic;
	interpreter* the_interpreter;
	interpreter* the_defaults_interpreter;
	
	char *argbuf;
	
public:
	
	Configuration( void )
	{
		ic = NULL;
		the_interpreter = NULL;
		the_defaults_interpreter = NULL;
		
		argbuf = new char[1024];
	}
	
	void Initialize( std::string interpreter_name, std::string input_file, char defaults_string[] )
	{
		if( the_interpreter != NULL ) delete the_interpreter;
		if( the_defaults_interpreter != NULL ) delete the_defaults_interpreter;
		
		ic = get_interpreter()[interpreter_name];
		the_interpreter = ic->create( input_file );
		the_defaults_interpreter = ic->create( std::string(defaults_string), true);
	}

	~Configuration()
	{
		delete the_interpreter;
		delete the_defaults_interpreter;
		delete[] argbuf;
	}

	template< typename T >
	void GetScalar( T& val, const char* key, ... ) const
	{
		va_list argptr;
		va_start(argptr, key);
		va_end(argptr);
		vsprintf(argbuf,key,argptr);		

		std::string strval;
		int status = 0;
		status = the_interpreter->query(std::string(argbuf),strval);
		if( status != 1 ) {
		  
		  fprintf(stderr, "Did not find <%s> in parameter file, trying defaults.\n",argbuf);

		  status = the_defaults_interpreter->query(std::string(argbuf),strval);
		  if( status != 1 ) {
		    fprintf(stderr, "Did not find <%s> in defaults either.\n",argbuf);
		    throw std::runtime_error("parameter not found!");
		  }

		}

		convert<std::string,T>(strval,val);
	}
	
	
	template< class T >
	void GetArray( T* val, const char* key, ... ) const
	{
		va_list argptr;
		va_start(argptr, key);
		va_end(argptr);
		vsprintf(argbuf,key,argptr);
		
		std::vector<std::string> s;
		int status = 0;
		status = the_interpreter->query_list(std::string(argbuf),s);
		if( status != 1 ) {
		  
		  fprintf(stderr, "Did not find <%s> in parameter file, trying defaults.\n",argbuf);
		  
		  status = the_defaults_interpreter->query_list(std::string(argbuf),s);
		  if( status != 1 ) {
		    fprintf(stderr, "Did not find <%s> in defaults either.\n",argbuf);
		    throw std::runtime_error("parameter not found!");
		  }
		}
		  
		for( size_t i=0; i<s.size(); ++i )
		  {
		    convert<std::string,T>(s[i],*val);
		    ++val;
		  }
	}
		
	size_t Size( const char* key, ... ) const
	{
		va_list argptr;
		va_start(argptr, key);
		va_end(argptr);
		vsprintf(argbuf,key,argptr);
		
		return the_interpreter->size(std::string(argbuf));
		
	}
	
};



template<>
void Configuration::GetScalar<char*>( char*& str, const char* key, ... ) const;

template<>
void Configuration::GetArray<char*>( char** val, const char* key, ... ) const;

/*********************************************************************/
/*********************************************************************/


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/


#endif
