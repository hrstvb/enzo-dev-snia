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

void EnableActiveParticleType(char *active_particle_type_name) {
    std::string dummy = std::string(active_particle_type_name);
    ActiveParticleType_info *my_type =
        get_active_particle_types()[dummy];
    if (my_type == NULL) {
        ENZO_FAIL("Unknown ParticleType");
    }
    EnabledActiveParticles[EnabledActiveParticlesCount++] = my_type;
    return;
}

