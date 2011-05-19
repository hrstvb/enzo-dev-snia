/***********************************************************************
/
/  AN EXAMPLE ACTIVE PARTICLE TYPE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/
************************************************************************/

#include <string.h>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
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
#include "EventHooks.h"
#include "ActiveParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_SampleParticle;
class SampleParticleBufferHandler;

class SampleParticleGrid : private grid {
  friend class ActiveParticleType_SampleParticle;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

class ActiveParticleType_SampleParticle : public ActiveParticleType
{
public:
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
};

int ActiveParticleType_SampleParticle::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  SampleParticleGrid *thisgrid =
    static_cast<SampleParticleGrid *>(thisgrid_orig);
  fprintf(stderr, "Checking formation of sample particles.\n");
  return 0;
}

void ActiveParticleType_SampleParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
}

class SampleParticleBufferHandler : public ParticleBufferHandler
{
  public:
    SampleParticleBufferHandler(int NumberOfParticles) { }
};

ParticleBufferHandler *ActiveParticleType_SampleParticle::AllocateBuffers(int NumberOfParticles)
{
    SampleParticleBufferHandler *handler = new SampleParticleBufferHandler(NumberOfParticles);
    return handler;
}


namespace {
  ActiveParticleType_info *SampleInfo = new ActiveParticleType_info
    ("SampleParticle", (&ActiveParticleType_SampleParticle::EvaluateFormation),
     (&ActiveParticleType_SampleParticle::DescribeSupplementalData),
     (&ActiveParticleType_SampleParticle::AllocateBuffers));
}
