/***********************************************************************
/
/  READS RADIATIVE TRANSFER PARAMETERS FROM INPUT FILE
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE: 
/
/  NOTE: modeled after CosmologyReadParameters
/
************************************************************************/

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* Set default parameter values. */

const char config_radiative_transfer_defaults[] = 
"### RADIATIVE TRANSFER INITIALIZATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    RadiativeTransfer: {\n"
"        RadiationPressure          = False;\n"
"        RadiationPressureScale     = 1.0;\n"
"        PhotonTime                 = 0;\n"
"        dtPhoton                   = -9999.0;\n"
"        SourceRadius               = 0;\n"
"        PropagationSpeedFraction   = 1.0;\n"
"        PropagationDistance        = 0.1;\n"
"        CoupledRateSolver          = True;\n"
"        OpticallyThinH2            = True;\n"
"        FluxBackgroundLimit        = 0.01;\n"
"        SplitPhotonRadius          = -99999.0;  # kpc\n"
"        RaysPerCell                = 5.1;\n"
"        InitialHEALPixLevel        = 3;\n"
"        PhotonEscapeRadius         = 0.0;  # kpc\n"
"        InterpolateField           = False;\n"
"        SourceClustering           = False;\n"
"        PhotonMergeRadius          = 10.0;\n"
"        TimestepVelocityLimit      = 100.0;  # km/s\n"
"        PeriodicBoundary           = False;\n"
"        FLDCallOnLevel             = 0;\n"
"        HIIRestrictedTimestep      = False;\n"
"        AdaptiveTimestep           = False;\n"
"        HydrogenOnly               = False;\n"
"        TraceSpectrum              = False;\n"
"        TraceSpectrumTable         = \"spectrum_table.dat\";\n"
"        SourceBeamAngle            = 30.0;\n"
"        LoadBalance                = False;\n"
"        OpticallyThinH2            = False;\n"
"    }\n"
"}\n";  


int RadiativeTransferReadParameters()
{
  int i;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  /* Set some defaults for non-parameter variables */
  for (i = 0; i < 4; i++) {
    EscapedPhotonCount[i] = 0.0;
    TotalEscapedPhotonCount[i] = 0.0;
  }
  PhotonEscapeFilename   = NULL;
  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;
  SourceClusteringTree = NULL;
  OldSourceClusteringTree = NULL;

  // This is how it should look eventually.
  //Param.UpdateDefaults(config_radiative_transfer_defaults);

  /* read parameters */
    
  Param.GetScalar(RadiationPressure, "Problem.RadiativeTransfer.RadiationPressure");
  Param.GetScalar(RadiationPressureScale, "Problem.RadiativeTransfer.RadiationPressureScale");
  Param.GetScalar(RadiativeTransferSourceRadius, "Problem.RadiativeTransfer.SourceRadius");
  Param.GetScalar(RadiativeTransferPropagationSpeedFraction, "Problem.RadiativeTransfer.PropagationSpeedFraction");
  Param.GetScalar(RadiativeTransferPropagationDistance, "Problem.RadiativeTransfer.PropagationDistance");
  Param.GetScalar(RadiativeTransferCoupledRateSolver, "Problem.RadiativeTransfer.CoupledRateSolver");
  Param.GetScalar(RadiativeTransferOpticallyThinH2, "Problem.RadiativeTransfer.OpticallyThinH2");
  Param.GetScalar(RadiativeTransferPeriodicBoundary, "Problem.RadiativeTransfer.PeriodicBoundary");
  Param.GetScalar(RadiativeTransferSplitPhotonRadius, "Problem.RadiativeTransfer.SplitPhotonRadius");
  Param.GetScalar(RadiativeTransferFluxBackgroundLimit, "Problem.RadiativeTransfer.FluxBackgroundLimit");
  Param.GetScalar(RadiativeTransferRaysPerCell, "Problem.RadiativeTransfer.RaysPerCell");
  Param.GetScalar(RadiativeTransferTimestepVelocityLimit, "Problem.RadiativeTransfer.TimestepVelocityLimit");
  Param.GetScalar(RadiativeTransferInitialHEALPixLevel, "Problem.RadiativeTransfer.InitialHEALPixLevel");
  Param.GetScalar(RadiativeTransferPhotonEscapeRadius, "Problem.RadiativeTransfer.PhotonEscapeRadius");
  Param.GetScalar(RadiativeTransferInterpolateField, "Problem.RadiativeTransfer.InterpolateField");
  Param.GetScalar(RadiativeTransferSourceClustering, "Problem.RadiativeTransfer.SourceClustering");
  Param.GetScalar(RadiativeTransferPhotonMergeRadius, "Problem.RadiativeTransfer.PhotonMergeRadius");
  Param.GetScalar(RadiativeTransferFLDCallOnLevel, "Problem.RadiativeTransfer.FLDCallOnLevel");
  Param.GetScalar(RadiativeTransferSourceBeamAngle, "Problem.RadiativeTransfer.SourceBeamAngle");
  Param.GetScalar(RadiativeTransferHIIRestrictedTimestep, "Problem.RadiativeTransfer.HIIRestrictedTimestep");
  Param.GetScalar(RadiativeTransferAdaptiveTimestep, "Problem.RadiativeTransfer.AdaptiveTimestep");
  Param.GetScalar(RadiativeTransferHydrogenOnly, "Problem.RadiativeTransfer.HydrogenOnly");
  Param.GetScalar(RadiativeTransferTraceSpectrum, "Problem.RadiativeTransfer.TraceSpectrum");
  Param.GetScalar(RadiativeTransferLoadBalance, "Problem.RadiativeTransfer.LoadBalance");
  Param.GetScalar(dtPhoton, "Problem.RadiativeTransfer.dtPhoton");

  Param.GetScalar(dummy, "Problem.RadiativeTransfer.TraceSpectrumTable");
  if( strlen(dummy) > 0 ) {
    RadiativeTransferTraceSpectrumTable = new char[MAX_LINE_LENGTH];
    strcpy(RadiativeTransferTraceSpectrumTable, dummy);
  }
  

  /* Error check */

  /* Check if H2 cooling is turned on for Lyman-Werner radiation. */

  if (RadiativeTransferOpticallyThinH2 && MultiSpecies < 2) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: optically thin Lyman-Werner radiation turned on "
	      "without H2 cooling.  Setting LW radiation OFF.\n");
    RadiativeTransferOpticallyThinH2 = FALSE;
  }

  /* Check if we're simultaneously using FLD and 1/r^2 Lyman-Werner radiation */

  if (RadiativeTransferOpticallyThinH2 && RadiativeTransferFLD) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: optically thin Lyman-Werner radiation and FLD "
	      "turned on.  Turning the optically thin radiation OFF.\n");
    RadiativeTransferOpticallyThinH2 = FALSE;
  }

  if (RadiativeTransferFLDCallOnLevel < 0) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: RadiativeTransferFLDCallOnLevel = %"ISYM
	      " cannot be negative!  Setting to 0.\n", RadiativeTransferFLDCallOnLevel);
    RadiativeTransferFLDCallOnLevel = 0;
  }


  // If RadiativeTransferFLD > 1, turn off RadiativeTransfer
  if (RadiativeTransferFLD > 1  &&  RadiativeTransfer) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: RadiativeTransferFLD > 1 cannot be used with "
	      "RadiativeTransfer.  Turning ray-tracing solver OFF.\n");
    RadiativeTransfer = FALSE;
  }


  // If RadiativeTransferFLD > 1, turn off RadiativeCooling
  if (RadiativeTransferFLD > 1  &&  RadiativeCooling) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: RadiativeTransferFLD > 1 cannot be used with "
	      "RadiativeCooling solver.  Turning RadiativeCooling OFF.\n");
    RadiativeCooling = 0;
  }

  // If RadiativeTransferFLD > 1, reset RadiationFieldType (if necessary)
  if (RadiativeTransferFLD > 1  &&  (RadiationFieldType != 0)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning: RadiativeTransferFLD > 1 cannot be used with "
	      "RadiationFieldType != 0.  Resetting RadiationFieldType.\n");
    RadiationFieldType = 0;
  }

  // If RadiativeTransferFLD > 1, ensure that ImplicitProblem > 0
  if (RadiativeTransferFLD > 1  &&  (ImplicitProblem == 0)) 
    ENZO_FAIL("Error: RadiativeTransferFLD > 1 requires ImplicitProblem > 0!")

  delete [] dummy;

  return SUCCESS;
}
