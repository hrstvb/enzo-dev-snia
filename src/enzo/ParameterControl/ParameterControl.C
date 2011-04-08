#include "ParameterControl.h"



template<>
void Configuration::GetScalar<char*>( char*& str, const char* key, ... ) const
{
	va_list argptr;
	va_start(argptr, key);
	va_end(argptr);
	vsprintf(argbuf,key,argptr);
	
	std::string strval;
	the_interpreter->query(std::string(argbuf),strval);
	
	strcpy( str, strval.c_str() );
}



template<>
void Configuration::GetArray<char*>( char** val, const char* key, ... ) const
{
	va_list argptr;
	va_start(argptr, key);
	va_end(argptr);
	vsprintf(argbuf,key,argptr);
	
	std::vector<std::string> s;
	
	the_interpreter->query_list(std::string(argbuf),s);
	
	for( size_t i=0; i<s.size(); ++i )
	{
		strcpy(*val,s[i].c_str());
		++val;
	}
}



std::map< std::string, interpreter_creator*>& get_interpreter()
{
	static std::map< std::string, interpreter_creator* > interpreter_map;
	return interpreter_map;
}


