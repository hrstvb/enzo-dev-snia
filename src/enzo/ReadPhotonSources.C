#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#include <stdlib.h>
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
#include "LevelHierarchy.h"

#define RT_ENERGY_BINS 4
#define MAX_RT_ENERGY_BINS 100

/* Set default parameter values. */

const char config_radiative_transfer_photon_sources_defaults[] = 
"### RADIATIVE TRANSFER PHOTON SOURCES DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    RadiativeTransfer: {\n"
"        PhotonSources = [\"Source1\"];\n"
"        Source1: {\n"
"            Type = 1;  # isotropic\n"
"            Luminosity = 0.0;\n"
"            LifeTime = 0.;"
"            CreationTime = -99999.0;\n"
"            RampTime = 0.0;\n"
"            SED = [1.0];\n"
"            Energy = [14.6];\n"
"            Position = [-99999.0, -99999.0, -99999.0];\n"
"        }\n"
"    }\n"
"}\n";  

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);

int ReadPhotonSources(FLOAT CurrentTime)
{

  int i, source, dim;

  int   PhotonTestNumberOfSources;
  int   PhotonTestSourceType[MAX_SOURCES];
  int   PhotonTestSourceEnergyBins[MAX_SOURCES];
  double PhotonTestSourceLuminosity[MAX_SOURCES];
  FLOAT PhotonTestSourcePosition[MAX_SOURCES][MAX_DIMENSION];
  float PhotonTestSourceLifeTime[MAX_SOURCES];
  float PhotonTestSourceCreationTime[MAX_SOURCES];
  float PhotonTestSourceRampTime[MAX_SOURCES];
  float *PhotonTestSourceSED[MAX_SOURCES];
  float *PhotonTestSourceEnergy[MAX_SOURCES];

  // Update the parameter config to include the local defaults. Note
  // that this does not overwrite values previously specified.
  Param.Update(config_radiative_transfer_photon_sources_defaults);

  // We need to update the Source1.CreationTime and Position
  // parameters in the defaults settings, since they depends on
  // runtime parameters (CreationTime, DomainLeftEdge,
  // DomainRightEdge).

  Param.SetScalar(CurrentTime, "Problem.RadiativeTransfer.Source1.CreationTime");
  FLOAT DefaultSourcePosition[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    DefaultSourcePosition[dim] = 0.5*(DomainLeftEdge[dim] +
				      DomainRightEdge[dim]);
  Param.SetArray(DefaultSourcePosition, MAX_DIMENSION, "Problem.RadiativeTransfer.Source1.Position");


  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, CurrentTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  char *numbers;
  char *delims = (char*) " ";
  char *value;
  int count;
  bool EnergyBinsDefined = false;


  /* read parameters */
  PhotonTestNumberOfSources = Param.Size("Problem.RadiativeTransfer.PhotonSources");
  if (PhotonTestNumberOfSources > MAX_SOURCES) {
    ENZO_VFAIL("You've exceeded the maximum number of RadiativeTransfer.PhotonSources (%d)!\n",MAX_SOURCES)
  }

  char PhotonSourceNames[MAX_LINE_LENGTH][MAX_SOURCES];
  Param.GetArray(PhotonSourceNames, "Problem.RadiativeTransfer.PhotonSources");

  for (i = 0; i < PhotonTestNumberOfSources; i++) {
    if (debug)
      fprintf(stdout, "ReadPhotonSources: Reading Parameters of Source %"ISYM"...\n", i);

    Param.GetScalar(PhotonTestSourceType[i], "Problem.RadiativeTransfer.%s.Type",PhotonSourceNames[i]);
    Param.GetArray(PhotonTestSourcePosition[i], "Problem.RadiativeTransfer.%s.Position",PhotonSourceNames[i]);
    Param.GetScalar(PhotonTestSourceLuminosity[i], "Problem.RadiativeTransfer.%s.Luminosity",PhotonSourceNames[i]);
    Param.GetScalar(PhotonTestSourceCreationTime[i], "Problem.RadiativeTransfer.%s.CreationTime",PhotonSourceNames[i]);
    Param.GetScalar(PhotonTestSourceLifeTime[i], "Problem.RadiativeTransfer.%s.LifeTime",PhotonSourceNames[i]);
    Param.GetScalar(PhotonTestSourceRampTime[i], "Problem.RadiativeTransfer.%s.RampTime",PhotonSourceNames[i]);
    
    PhotonTestSourceSED[i] = new float[MAX_RT_ENERGY_BINS];
    PhotonTestSourceEnergy[i] = new float[MAX_RT_ENERGY_BINS];

    PhotonTestSourceEnergyBins[i] = Param.Size("Problem.RadiativeTransfer.%s.SED",PhotonSourceNames[i]);
    Param.GetArray(PhotonTestSourceSED[i], "Problem.RadiativeTransfer.%s.SED",PhotonSourceNames[i]);
    Param.GetArray(PhotonTestSourceEnergy[i], "Problem.RadiativeTransfer.%s.Energy",PhotonSourceNames[i]);

  }


  // normalize SED

  float totSED;
  for (source = 0; source < PhotonTestNumberOfSources; source++) {
    totSED = 0.;  
    for (i=0; i<PhotonTestSourceEnergyBins[source]; i++) 
      totSED += PhotonTestSourceSED[source][i];
    for (i=0; i<PhotonTestSourceEnergyBins[source]; i++) 
      PhotonTestSourceSED[source][i] /= totSED;
  }

  // list head:
  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;
  for (i=0; i<PhotonTestNumberOfSources; i++) {
    if (debug) fprintf(stdout, "ReadPhotonSources: %"ISYM" %"GSYM" %"GSYM" %"GSYM"\n", 
		       i, PhotonTestSourceLuminosity[i], TimeUnits, LengthUnits);
    PhotonTestSourceLuminosity[i] *= TimeUnits/pow(LengthUnits,3);
    if (debug) fprintf(stdout, "ReadPhotonSources: %"ISYM"  %"GSYM"\n", 
		       i, PhotonTestSourceLuminosity[i]);
    RadiationSourceEntry *RadSources;
    RadSources = new RadiationSourceEntry;
    RadSources->PreviousSource = GlobalRadiationSources;
    RadSources->NextSource     = GlobalRadiationSources->NextSource;
    RadSources->SuperSource    = NULL;
    RadSources->GridID         = INT_UNDEFINED;
    RadSources->GridLevel      = INT_UNDEFINED;
    RadSources->Type           = PhotonTestSourceType[i]; 
    RadSources->Luminosity     = PhotonTestSourceLuminosity[i]; 
    RadSources->LifeTime       = PhotonTestSourceLifeTime[i]; 
    RadSources->RampTime       = PhotonTestSourceRampTime[i]; 
    RadSources->CreationTime   = PhotonTestSourceCreationTime[i]; 
    RadSources->Position       = new FLOAT[3];
    RadSources->Position[0]    = PhotonTestSourcePosition[i][0]; 
    RadSources->Position[1]    = PhotonTestSourcePosition[i][1]; 
    RadSources->Position[2]    = PhotonTestSourcePosition[i][2]; 
    RadSources->EnergyBins     = PhotonTestSourceEnergyBins[i];
    RadSources->Energy         = new float[PhotonTestSourceEnergyBins[i]];
    RadSources->SED            = new float[PhotonTestSourceEnergyBins[i]];
    RadSources->AddedEmissivity = false;
    for (int j=0; j<PhotonTestSourceEnergyBins[i]; j++){
      RadSources->Energy[j] = PhotonTestSourceEnergy[i][j];
      RadSources->SED[j]    = PhotonTestSourceSED[i][j];
    }

    if (RadSources->Type != Isotropic && RadSources->Type != Beamed &&
	RadSources->Type != Episodic) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "PhotonTestSourceType must be 1, -2, -3.\n",
		"\tChanging to 1 (isotropic)\n");
      RadSources->Type = Isotropic;
    }

    GlobalRadiationSources->NextSource = RadSources;

  }  

  /* Delete allocated memory for temporary (for I/O) photon sources */

  for (source = 0; source < PhotonTestNumberOfSources; source++) {
    delete [] PhotonTestSourceEnergy[source];
    delete [] PhotonTestSourceSED[source];
  }

  /* Create tree that clusters the sources if requested */

  /* While creating tree (type SuperSource), compute position of the
     super source in each leaf. */
  
  if (RadiativeTransferSourceClustering == TRUE)
    if (CreateSourceClusteringTree(NULL, NULL, NULL) == FAIL) {
      ENZO_FAIL("Error in CreateSourceClusteringTree.\n");

    }
  
  return SUCCESS;

}
