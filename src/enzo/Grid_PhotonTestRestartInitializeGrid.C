/***********************************************************************
/
/  GRID CLASS (RE-INITIALIZE FOR PUTTING IN SINK PARTICLES)
/
/  written by: Elizabeth Harper-Clark
/  date:       March, 2010
/  modified1:
/
/  PURPOSE: based on SupernovaeRestartInitialize
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::PhotonTestRestartInitialize(int level, int *NumberOfCellsSet)
{

  if (RadiativeTransfer && (MultiSpecies < 1)) {
    ENZO_FAIL("Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
  }

  NumberOfPhotonPackages = 0;
  PhotonPackages-> NextPackage= NULL;

  /* Initialize radiation fields - not needed in restart?? */

//   if (this->InitializeRadiativeTransferFields() == FAIL) {
//     ENZO_FAIL("\nError in InitializeRadiativeTransferFields.\n");
//   }



  return SUCCESS;
}
