/***********************************************************************
/
/  EVOLVE PHOTONS FUNCTION
/
/  written by: Tom Abel
/  date:       May 2004
/  modified1:  November 2005 by John Wise (parallelized it)
/                
/
/  PURPOSE:
/    This routine is the main photon evolution function. 
/    From here we call first the routines that control the emission:
/    grid::Shine
/  while (stil_photons_to_update)
/    transport all photon packages local to a grid
/    communicate all surviving photon packages to their new parent grids. 
/  endwhile  
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"

/* function prototypes */
void my_exit(int status);
int CommunicationTransferPhotons(LevelHierarchyEntry *LevelArray[], 
				 ListOfPhotonsToMove **AllPhotons,
				 int &keep_transporting);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int InitiateKeepTransportingCheck(int keep_transporting);
int StopKeepTransportingCheck();
int InitializePhotonCommunication();
int KeepTransportingCheck(int &keep_transporting);
RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);
void PrintSourceClusteringTree(SuperSourceEntry *leaf);
int OutputEscapeFraction(void);
int CommunicationSyncNumberOfPhotons(LevelHierarchyEntry *LevelArray[]);
int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level);

/* EvolvePhotons function */
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, FLOAT GridTime, int level, int LoopTime)
{

  bool FirstTime = true;

  if (!RadiativeTransfer)
    return SUCCESS;

  /* Only call on the finest level */

  if (LevelArray[level+1] != NULL)
    return SUCCESS;

//  printf("GridTime = %f, PhotonTime = %f, dtPhoton = %g (Loop = %d)\n",
//	 GridTime, PhotonTime, dtPhoton, (GridTime >= PhotonTime));

  //if (dtPhoton < 0)
  //  return SUCCESS;

  //while (GridTime >= PhotonTime) {
  while (GridTime > PhotonTime) {

    /* Recalculate timestep if this isn't the first loop.  We already
       did this in RadiativeTransferPrepare */

    FLOAT dtLevelAbove;
    if (!FirstTime) {
      dtLevelAbove = LevelArray[level]->GridData->ReturnTimeStep();
      RadiativeTransferComputeTimestep(LevelArray, MetaData, dtLevelAbove, level);
    }

    if (debug && LoopTime == TRUE)
      printf("EvolvePhotons[%"ISYM"]: dt = %"GSYM", Time = %"FSYM", ", 
	     level, dtPhoton, PhotonTime);
      
    /* Declarations */

    grid *Helper;
    int RefinementFactors[MAX_DIMENSION];

    /* Create an array (Grids) of all the grids. */

    typedef HierarchyEntry* HierarchyEntryPointer;
    HierarchyEntry **Grids;
    HierarchyEntry **Parents;
    LevelHierarchyEntry *Temp;
    RadiationSourceEntry *RS;
    int GridNum = 0, value, i, proc, lvl;
    int NumberOfGrids = 0;  
    int NumberOfSources, NumberOfSourcesInLoop, SourcesCompleted;
    int SourcesCompletedInLoop;

    /* delete source if we are passed (or before) their lifetime */

    RS = GlobalRadiationSources->NextSource;
    NumberOfSources = 0;
    SourcesCompleted = 0;

    while (RS != NULL) {
      if ( (RS->CreationTime + RS->LifeTime) < PhotonTime ||
	   (RS->CreationTime > PhotonTime + dtPhoton) ) {
	if (debug) {
	  fprintf(stdout, "\nEvolvePhotons: Deleted Source on lifetime limit \n");
	  fprintf(stdout, "EvolvePhotons:  %"GSYM" %"GSYM" %"GSYM" \n",
		  RS->CreationTime, RS->LifeTime, PhotonTime);
	}
	RS = DeleteRadiationSource(RS);
      } else {
	NumberOfSources++;                 // count sources
	RS = RS->NextSource;
      }
    }

    if (debug) fprintf(stdout, "%"ISYM" SRC(s)\n", NumberOfSources);

    int Rank, Dims[MAX_DIMENSION];
    FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  
    /* Initialize radiation fields */  

    for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->InitializeRadiativeTransferFields();

    for (i = 0; i < 4; i++)
      EscapedPhotonCount[i] = 0.0;

    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
      FieldsToInterpolate[i] = FALSE;

    if (NumberOfSources == 0) {
      PhotonTime += dtPhoton;
      continue;
    }    

    /* Create tree that clusters the sources if requested.  While
       creating tree (type SuperSource), compute position of the super
       source in each leaf. */

    if (RadiativeTransferSourceClustering == TRUE)
      CreateSourceClusteringTree(NULL, NULL, LevelArray);

    /* We need a list of root grids for transport between grids */

    HierarchyEntry **Temp0;
    int nGrids0 = GenerateGridArray(LevelArray, 0, &Temp0);
    grid **Grids0 = new grid*[nGrids0];
    for (i = 0; i < nGrids0; i++)
      Grids0[i] = Temp0[i]->GridData;

    RS = GlobalRadiationSources->NextSource;

    /* Because the ray marking field only holds N bits of information,
       we loop over that number of sources. */

    NumberOfSourcesInLoop = 8*sizeof(int);

    while (SourcesCompleted < NumberOfSources) {

    SourcesCompletedInLoop = 0;
 
    while (RS != NULL && SourcesCompletedInLoop < NumberOfSourcesInLoop) {
      int Continue = 1;
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; (lvl >= 0 && Continue); lvl--) {
	for (Temp = LevelArray[lvl]; (Temp && Continue); 
	     Temp = Temp->NextGridThisLevel)
	  if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())

	    // Find grid
	    if (Temp->GridData->PointInGrid(RS->Position)) {
	      Temp->GridData->Shine(RS);
	      Continue = FALSE; // do not continue with this source
	    } // If source in grid

	/* For MPI, communicate to minimum value of Continue to ensure
	   that the source's host grid was found. */

	Continue = CommunicationMinValue(Continue);

	if (lvl == 0 && Continue) {  // this should never happen ... 
	  fprintf(stderr, "Could not find grid for source %x: "
		  "Pos: %"FSYM" %"FSYM" %"FSYM"\n",
		  RS, RS->Position[0], RS->Position[1], RS->Position[2]);
	  ENZO_FAIL("");
	}
      }    // Loop through levels 
      RS = RS->NextSource;
      SourcesCompletedInLoop++;
    }    // ENDWHILE still sources 

    SourcesCompleted += SourcesCompletedInLoop;

#ifdef USE_MPI
    if (RadiativeTransferInterpolateField)
      CommunicationAllReduceValues(FieldsToInterpolate,
				   MAX_NUMBER_OF_BARYON_FIELDS, MPI_MAX);
#endif /* USE_MPI */  

    /* Initialize ray marker */

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->InitializeRayMarker();

    /* Evolve all photons by fixed timestep. */
  
    ListOfPhotonsToMove *PhotonsToMove = new ListOfPhotonsToMove;
    PhotonsToMove->NextPackageToMove = NULL;

    int keep_transporting = 1;
    int ThisProcessor;

    /* Initialize nonblocking MPI routine */

    InitializePhotonCommunication();

    /* Transport the rays! */

    while (keep_transporting) {
      keep_transporting = 0;
      PhotonsToMove->NextPackageToMove = NULL;
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--) {

	NumberOfGrids = GenerateGridArray(LevelArray, lvl, &Grids);
	for (GridNum = 0; GridNum < NumberOfGrids; GridNum++) {

	  if (Grids[GridNum]->ParentGrid != NULL) 
	    Helper = Grids[GridNum]->ParentGrid->GridData;
	  else
	    Helper = NULL;

	  Grids[GridNum]->GridData->TransportPhotonPackages
	    (lvl, &PhotonsToMove, GridNum, Grids0, nGrids0, Helper, 
	     Grids[GridNum]->GridData);

	} // ENDFOR grids

	delete [] Grids;

      }                          // loop over levels

      if (PhotonsToMove->NextPackageToMove != NULL)
	keep_transporting = 1;

      /* Check if there are any photons leaving this grid.  If so,
	 move them. */
      
      CommunicationTransferPhotons(LevelArray, &PhotonsToMove, 
				   keep_transporting);

      /* Receive keep_transporting messages and take the MAX */

#ifdef NONBLOCKING
      InitiateKeepTransportingCheck(keep_transporting);
      KeepTransportingCheck(keep_transporting);
#else /* NON_BLOCKING */
      keep_transporting = CommunicationMaxValue(keep_transporting);
#endif

    }                           //  end while keep_transporting

    delete PhotonsToMove;

    } // ENDWHILE (SourcesCompleted < NumberOfSources)

    //  StopKeepTransportingCheck();

    /* Move all finished photon packages back to their original place,
       PhotonPackages */

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->MoveFinishedPhotonsBack();

    OutputEscapeFraction();

    PhotonTime += dtPhoton;

    /* Delete ray marker */

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->DeleteRayMarker();

    delete [] Grids0;
    delete [] Temp0;

    /* Coupled rate & energy solver */

    float dtMin, dtGrid;
    float dtLastLevel = 1e20;

    int debug_store = debug;
    debug = FALSE;

    /* Divide the photo-ionization and photo-heating rates by the
       number of particles (rho * dx^3) */

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE)
	  Temp->GridData->FinalizeRadiationFields();

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE) {

	  if (RadiativeTransferCoupledRateSolver && 
	      RadiativeTransferOpticallyThinH2)
	    Temp->GridData->AddH2Dissociation(AllStars);

	  if (RadiativeTransferCoupledRateSolver)
	    Temp->GridData->SolveCoupledRateEquations();

	  if (RadiativeTransferCoupledRateSolver &&
	      RadiativeTransferInterpolateField)
	    Temp->GridData->DeleteInterpolatedFields();

	} /* ENDIF radiation */

    /* For the non-coupled (i.e. cells without radiation) rate & energy
       solver, we have to set the H2 dissociation rates */

    if (RadiativeTransferOpticallyThinH2)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->RadiationPresent() == FALSE)
	    Temp->GridData->AddH2Dissociation(AllStars);

    debug = debug_store;

    /* We don't rely on the count NumberOfPhotonPackages here, so they
       aren't synchronized across processors.  But in RebuildHierarchy, 
       this number is needed.  Synchronize them now. */

    CommunicationSyncNumberOfPhotons(LevelArray);

    /* If we're using the HII restricted timestep, get the global
       maximum kph in I-fronts. */

    if (RadiativeTransferHIIRestrictedTimestep) {
      float LocalMaximumkph = -1e20;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  LocalMaximumkph = max(LocalMaximumkph,
				Temp->GridData->ReturnMaximumkphIfront());
      LocalMaximumkph = CommunicationMaxValue(LocalMaximumkph);
      MetaData->GlobalMaximumkphIfront = LocalMaximumkph;
    }

#ifdef DEBUG
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->ErrorCheckPhotonNumber(lvl);
#endif

    if (!LoopTime)
      break;
    
    FirstTime = false;

  } // ENDWHILE GridTime >= PhotonTime

  return SUCCESS;

}
