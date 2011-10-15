/***********************************************************************
/
/  KRAVTSOV (2003) STAR FORMATION
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
#include "phys_constants.h"

#ifdef NEW_CONFIG

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

/* Set default parameter values. */

const char config_kravtsov_particle_defaults[] = 
"### KRAVTSOV STAR PARTICLE DEFAULTS ###\n"
"\n"
"Physics: {\n"
"    ActiveParticles: {\n"
"        Kravtsov: {\n"
"            DensityThreshold = 1e6; # [particles per proper cm^3]\n"
"            StarFormationTimeConstant = 4.0e9; # [years]\n"
"            MinimumStarMass = 1.0e9; # [Msun]\n"
"        };\n"
"    };\n"
"};\n";

#endif


/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class ActiveParticleType_Kravtsov;
class KravtsovParticleBufferHandler;

class KravtsovGrid : private grid {
    friend class ActiveParticleType_Kravtsov;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

class ActiveParticleType_Kravtsov : public ActiveParticleType
{
public:
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static ParticleBufferHandler *AllocateBuffers(int NumberOfParticles);
  static int InitializeParticleType();
  static int EvaluateFeedback(grid *thisgrid_orig);

  static float DensityThreshold, StarFormationTimeConstant, MinimumStarMass;

private:
  float Metallicity;

};

float ActiveParticleType_Kravtsov::DensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::StarFormationTimeConstant = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::MinimumStarMass = FLOAT_UNDEFINED;

int ActiveParticleType_Kravtsov::InitializeParticleType() {

#ifdef NEW_CONFIG

  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite any values previously specified.
  Param.Update(config_kravtsov_particle_defaults);

  // Retrieve parameters from Param structure
  Param.GetScalar(DensityThreshold, "Physics.ActiveParticles.Kravtsov.DensityThreshold");
  Param.GetScalar(StarFormationTimeConstant, "Physics.ActiveParticles.Kravtsov.StarFormationTimeConstant");
  Param.GetScalar(MinimumStarMass, "Physics.ActiveParticles.Kravtsov.MinimumStarMass");

#else

  DensityThreshold = 0.0;
  StarFormationTimeConstant = 0.0;
  MinimumStarMass = 0.0;
  
#endif

  return SUCCESS;
}


int ActiveParticleType_Kravtsov::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  KravtsovGrid *tg =
    static_cast<KravtsovGrid *>(thisgrid_orig);


  /* Make it pretty */

  float *density = tg->BaryonField[supp_data.DensNum];
  float *velx = tg->BaryonField[supp_data.Vel1Num];
  float *vely = tg->BaryonField[supp_data.Vel2Num];
  float *velz = tg->BaryonField[supp_data.Vel3Num];

  float CellWidthTemp = float(tg->CellWidth[0][0]);

  bool HasMetalField = (supp_data.MetalNum != -1 ||
			supp_data.ColourNum != -1);

  int GridDimension[3] = {tg->GridDimension[0],
                          tg->GridDimension[1],
                          tg->GridDimension[2]};

  float gasfrac, starmass, densthresh, timeconstant;
  int i,j,k,index;
  int NumberOfNewParticles = 0;


  // calculate density threshold.  odthresh is in proper particles per
  // cc and (d1/mproton) gives the mean density of the universe in
  // particles/cm^3 (assuming hydrogen is dominant)

  densthresh = DensityThreshold / (supp_data.DensityUnits / mh);


  // calculate time constant for star formation.  This assumes that
  // the user input is in units of years
  
  timeconstant = StarFormationTimeConstant * 3.156e7 / supp_data.TimeUnits;


  // for each zone, : "star" particle is created if the density
  // exceeds some threshold and this is the highest level of
  // refinement.  That's it.

  for (k = tg->GridStartIndex[2]; k <= tg->GridEndIndex[2]; k++) {
    for (j = tg->GridStartIndex[1]; j <= tg->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(tg->GridStartIndex[0], j, k);
      for (i = tg->GridStartIndex[0]; i <= tg->GridEndIndex[0]; i++, index++) {
	
	// 0. If no more room for particles, quit.
	if (supp_data.NumberOfNewParticles >=
	    supp_data.MaxNumberOfNewParticles)
          continue;
	
	// 1. Finest level of refinement
	if (tg->BaryonField[tg->NumberOfBaryonFields][index] != 0.0) 
	  continue;
	
	// 2. Density greater than threshold
	if (density[index] < densthresh)
	  continue;
	
	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */
	
	ActiveParticleType_Kravtsov *np = new ActiveParticleType_Kravtsov();
	supp_data.NewParticles[supp_data.NumberOfNewParticles++] = np;

	// Make sure that we never give put than 90% of the cell's mass into a star particle
	gasfrac = min( 0.9, tg->dtFixed / timeconstant );
	
	// Calculate star mass in solar masses.  If this is less than
	// the user-defined threshold mass, do NOT make a star in this
	// cell.  This is not exactly in keeping with the spirit of
	// the Kravtsov algorithm, and is somewhat degenerate with the
	// density threshold, but we really don't want millions and
	// millions of star particles.

	starmass = gasfrac * density[index] * supp_data.DensityUnits * pow(supp_data.LengthUnits * CellWidthTemp, 3) /  SolarMass;

	// Do not allow stars with mass less than MinimumStarMass
	if (starmass < MinimumStarMass )
	  continue;

	np->Mass =  starmass / supp_data.MassUnits;

	np->type = Kravtsov;
	np->BirthTime = tg->Time;
	
	np->pos[0] = tg->CellLeftEdge[0][i] + 0.5*tg->CellWidth[0][i];
	np->pos[1] = tg->CellLeftEdge[1][j] + 0.5*tg->CellWidth[1][j];
	np->pos[2] = tg->CellLeftEdge[2][k] + 0.5*tg->CellWidth[2][k];

	/*
	  Star velocities averaged over multiple cells to avoid
	  "runaway star particle" phenomenon imethod = 2 is zeus,
	  otherwise PPM
	*/

	float *tvel = tg->AveragedVelocityAtCell(index, supp_data.DensNum,
						  supp_data.Vel1Num);
	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];

	/* Set the metallicity */

	if (HasMetalField)
	  np->Metallicity = supp_data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
  
  return NumberOfNewParticles;
}

// Pop III feedback
int ActiveParticleType_Kravtsov::EvaluateFeedback(grid *thisgrid_orig)
{
  return SUCCESS;
}

void ActiveParticleType_Kravtsov::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
  flags.MetalField = true;
}

class KravtsovParticleBufferHandler : public ParticleBufferHandler
{
  public:
    KravtsovParticleBufferHandler(int NumberOfParticles) { }
};

ParticleBufferHandler *ActiveParticleType_Kravtsov::AllocateBuffers(int NumberOfParticles)
{
    KravtsovParticleBufferHandler *handler = new KravtsovParticleBufferHandler(NumberOfParticles);
    return handler;
}

namespace {
    ActiveParticleType_info *SampleInfo = new ActiveParticleType_info(
            "Kravtsov", (&ActiveParticleType_Kravtsov::EvaluateFormation),
	    (&ActiveParticleType_Kravtsov::DescribeSupplementalData),
	    (&ActiveParticleType_Kravtsov::AllocateBuffers),
	    (&ActiveParticleType_Kravtsov::InitializeParticleType),
	    (&ActiveParticleType_Kravtsov::EvaluateFeedback) );

}
