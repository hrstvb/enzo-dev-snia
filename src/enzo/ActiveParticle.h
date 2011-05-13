/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels, Matt Turk
/  date:       May, 2011 (converting from Star to ActiveParticle)
/
/  PURPOSE:
/
************************************************************************/
#ifndef __ACTIVE_PARTICLE_H
#define __ACTIVE_PARTICLE_H

#include "typedefs.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarBuffer.h"

class ActiveParticleType
{
  public:

  protected:

  private:
  
};

/*!
 * @brief implements abstract factory design pattern for plug-ins
 */
struct ActiveParticleType_creator
{
    //! create an instance of a plug-in
    virtual ActiveParticleType * create( ) const = 0;
    
    //! destroy an instance of a plug-in
    virtual ~ActiveParticleType_creator() { }
};

typedef std::map<std::string, ActiveParticleType_creator *> ActiveParticleMap;

//! maps the name of a plug-in to a pointer of the factory pattern 
ActiveParticleMap &get_active_particle_types();

/*!
 * @brief concrete factory pattern for plug-ins
 */
template< class DerivedActiveParticleType >
struct ActiveParticle_creator_concrete : public ActiveParticleType_creator
{
    //! register the plug-in by its name
    ActiveParticle_creator_concrete( const std::string& active_particle_type_name )
    {
        get_active_particle_types()[ active_particle_type_name ] = this;
    }
    
    //! create an instance of the plug-in
    ActiveParticleType * create( ) const
    {
        return new DerivedActiveParticleType( );
    }
};

//! failsafe version to select the plug-in
ActiveParticleType *select_active_particle_type( std::string active_particle_type_name );


#endif

