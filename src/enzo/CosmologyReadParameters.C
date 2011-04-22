/***********************************************************************
/
/  READS COSMOLOGY PARAMETERS FROM INPUT FILE
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits);
 
int CosmologyReadParameters(FLOAT *StopTime, FLOAT *InitTime)
{
 
  int i, OutputNumber;
  FLOAT FinalRedshift, CurrentRedshift;
  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  /* Set defaults. */
 
  HubbleConstantNow    = 0.701;
  OmegaMatterNow       = 0.279;
  OmegaLambdaNow       = 0.721;
  ComovingBoxSize      = 64;
  MaxExpansionRate     = 0.01;
  InitialRedshift      = 20;
  FinalRedshift        = 0;
 
  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++) {
    CosmologyOutputRedshift[i]     = -1;  // Never!!
    CosmologyOutputRedshiftName[i] = NULL;
  }
 
  /* read parameters */
 
  

  Param.GetScalar(HubbleConstantNow,"Physics.Cosmology.HubbleConstantNow");
  Param.GetScalar(OmegaMatterNow, "Physics.Cosmology.OmegaMatterNow");
  Param.GetScalar(OmegaLambdaNow, "Physics.Cosmology.OmegaLambdaNow");
  Param.GetScalar(ComovingBoxSize, "Physics.Cosmology.ComovingBoxSize");
  Param.GetScalar(MaxExpansionRate, "Physics.Cosmology.MaxExpansionRate");
  Param.GetScalar(InitialRedshift, "Physics.Cosmology.InitialRedshift");
  Param.GetScalar(FinalRedshift, "Physics.Cosmology.FinalRedshift");
  Param.GetScalar(CurrentRedshift, "Internal.CosmologyCurrentRedshift");

  Param.GetArray(CosmologyOutputRedshift, "OutputControl.RedshiftDump.OutputRedshifts");

  int NumberOfOutputRedshiftNames = Param.Size("OutputControl.RedshiftDump.OutputRedshiftNames");
  char OutputRedshiftNames[MAX_LINE_LENGTH][MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
  Param.GetArray(OutputRedshiftNames, "OutputControl.RedshiftDump.OutputRedshiftNames");
  for (i = 0; i < NumberOfOutputRedshiftNames; i++) {
    CosmologyOutputRedshiftName[i] = new char[MAX_LINE_LENGTH];
    strcpy(CosmologyOutputRedshiftName[i], OutputRedshiftNames[i]);
  }

 
  /* Initialize by finding the time at the initial redshift. */
 
  if (CosmologyComputeTimeFromRedshift(InitialRedshift,
				       &InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in ComputeTimeFromRedshift.\n");
  }
  if (*InitTime == 0.0)
    *InitTime = InitialTimeInCodeUnits;
 
  /* Now find the time at the end of the simulation. */
 
  if (CosmologyComputeTimeFromRedshift(FinalRedshift, StopTime) == FAIL) {
    ENZO_FAIL("Error in ComputeTimeFromRedshift.\n");
  }
 
  /* Convert the output redshifts into time, for later convenience. */
 
  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
    if (CosmologyOutputRedshift[i] != -1)
      CosmologyComputeTimeFromRedshift(CosmologyOutputRedshift[i],
				       &CosmologyOutputRedshiftTime[i]);
 
  /* Convert the time action redshift into time. */
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++)
    if (TimeActionRedshift[i] != -1)

      CosmologyComputeTimeFromRedshift(TimeActionRedshift[i],
				       &TimeActionTime[i]);

  delete [] dummy;
 
  return SUCCESS;
}
