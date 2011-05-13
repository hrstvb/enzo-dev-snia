/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels,
/              Matt Turk
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
    /* Several pure virtual functions */

    /* This should return the number of new star particles created, and should
     * create them. */
    
  protected:

  private:
  
};

//! maps the name of a plug-in to a pointer of the factory pattern
class ActiveParticleType_info;
typedef std::map<std::string, ActiveParticleType_info *> ActiveParticleMap;

ActiveParticleMap &get_active_particle_types();

ActiveParticleType *select_active_particle_type( std::string active_particle_type_name );

class ActiveParticleType_info
{
    public:
       
       /* We will add more functions to this as necessary */
       ActiveParticleType_info(
           std::string this_name,
           int (*ffunc)(grid *thisgrid_orig) ){
        this->formation_function = ffunc;
        get_active_particle_types()[this_name] = this;
       }

       int (*formation_function)(grid *thisgrid_orig);
};

#endif

