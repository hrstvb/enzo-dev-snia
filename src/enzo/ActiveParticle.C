/***********************************************************************
/
/  PROBLEM TYPE CLASS
/
/  written by: Matthew Turk, Oliver Hahn
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#include "ActiveParticle.h"

ActiveParticleMap& get_active_particle_types()
{
    static ActiveParticleMap active_particle_type_map;
    return active_particle_type_map;
}

/* This takes a string, grabs the (static) plugin map defined above, and
   returns the plugin creator for that. */
ActiveParticleType *select_active_particle_type( std::string active_particle_type_name)
{
    ActiveParticleType_creator *ept_creator = get_active_particle_types()
            [ active_particle_type_name ];

    /* Simply throw an error if no such plugin exists... */

    if( !ept_creator )
    {   
        ActiveParticleMap mymap = get_active_particle_types();

        for (ActiveParticleMap::const_iterator it 
                = mymap.begin(); it != mymap.end(); ++it) {
          std::cout << "Available: " << it->first << std::endl;
        }
        ENZO_FAIL("Unknown Active Particle plug-in.");
    }

    ActiveParticleType *ptype = ept_creator->create();

    return ptype;

}
