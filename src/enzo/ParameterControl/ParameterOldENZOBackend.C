#include <map>


#include "ParameterBackends.h"
#include "SimpleParameterParser.h"

class enzo20_backend : public interpreter{
protected:
	
	//enzo21_backend new_backend_;
	config_file cf_;
	std::map<std::string,std::string> old_new_map_;
	
	void add( std::string oldkey, std::string newkey )
	{
		old_new_map_.insert( std::make_pair( newkey,oldkey ) );
	}
	
public:
	
	explicit enzo20_backend( std::string fname )
	: interpreter( fname ), cf_( fname )
	{
		
#define _BEGIN_MAP_ add(
#define _TO_ ,
#define _END_ ); add(
#define _END_MAP_ "empty","empty");
		
#include "ParameterMap.inc"
		
#undef _BEGIN_MAP_
#undef _TO_
#undef _END_
#undef _END_MAP_
	}
	
	void query( const std::string key, std::string& ret )
	{ 
		ret = cf_.getValue<std::string>( old_new_map_[key] );  
	}
	
	void query_list( const std::string key, std::vector< std::string >& ret )
	{
		std::string tmp = cf_.getValue<std::string>( old_new_map_[key] );
		cf_.evaluate_list<std::string>( tmp, std::back_inserter( ret ) );
	}
	
	size_t size( const std::string key )
	{
		std::vector<std::string> tmp;
		query_list( key, tmp );
		return tmp.size();
	}
	
};

namespace{
	interpreter_creator_concrete< enzo20_backend > c01("enzo2.0");
}
