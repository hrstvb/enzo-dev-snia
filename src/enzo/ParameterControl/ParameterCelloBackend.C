/* new cello style parameter file wrapper */

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include "ParameterBackends.h"

#include "parameters/parameters.hpp"

class enzo21_backend : public interpreter
{
protected:
	std::string fname_;
	Parameters * param_;
	
	
	void group_splitter( const std::string& instr, std::vector<std::string>& strvec )
	{
		std::string str( instr );
		size_t i;
		
		while( (i=str.find_first_of('.')) < str.size() )
		{
			strvec.push_back( str.substr(0,i) );
			str = str.substr(i+1);
		}
		
		strvec.push_back( str );
	}
	
	void open_group( const std::string& key, std::string& element_key )
	{
		if( key.find_first_of('.') == key.size() )
		{
			element_key = key;
			return;
		}
		
		std::vector< std::string > grp_list;
		group_splitter( key, grp_list );
		
		if( grp_list.size() > 3 ) 
		{	
			throw std::runtime_error("More than two level nesting currently not supported by backend!");
		}
		
		param_->set_current_group( grp_list[0] );
		
		if( grp_list.size() > 2 )
			param_->set_current_subgroup( grp_list[1] );
		
		element_key = grp_list[grp_list.size()-1];
		
	}
	
	/*void make_lower(  std::string& data )
	 {
	 std::transform(data.begin(), data.end(), data.begin(), ::tolower);
	 }*/
	
public:
	explicit enzo21_backend( std::string fname )
	: interpreter(fname), fname_(fname)
	{
		
		param_ = new Parameters();
		param_->read(fname_.c_str());
		
	}
	
	
	~enzo21_backend()
	{
		delete param_;
	}
	
	void query( const std::string key, std::string& ret )
	{   
		std::string element_key;
		open_group( key, element_key );
		ret = param_->raw_string( element_key.c_str() );
	}
	
	void query_list( const std::string key, std::vector< std::string >& ret )
	{
		std::string element_key;
		open_group( key, element_key );
		
		int nelem = param_->list_length( element_key );
		for( int i=0; i<nelem; ++i )
			ret.push_back( param_->list_raw_string( i, element_key.c_str() ) );
	}
	
	size_t size( const std::string key )
	{
		param_->set_current_group("List");
		std::string element_key;
		open_group( key, element_key );
		//make_lower( element_key );
		return (size_t) param_->list_length( element_key );
	}
};


namespace{
	interpreter_creator_concrete< enzo21_backend > c00("enzo2_cello");
}
	