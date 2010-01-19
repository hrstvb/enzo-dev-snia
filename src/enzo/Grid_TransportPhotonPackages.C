#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (TRANSPORT PHOTON PACKAGES)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This is the heart of the radiative transfer algorithm.
/    On each Grid we initialize photo and heating rates and then call
/    WalkPhotonPackage so all photon packages are transported along their
/    own directions and the photo-ionization and heating rates on 
/    on the grid are updated on the fly. 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#ifdef _OPENMP
#include "omp.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

#ifdef CONFIG_BFLOAT_4
#define ROUNDOFF 1e-6
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
PhotonPackageEntry *PopPhoton(PhotonPackageEntry * &Node);
PhotonPackageEntry *DeletePhotonPackage(PhotonPackageEntry *PP);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::TransportPhotonPackages(int level, ListOfPhotonsToMove **PhotonsToMove, 
				  int GridNum, grid **Grids0, int nGrids0, 
				  grid *ParentGrid, grid *CurrentGrid)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || MultiSpecies < 1) 
    return SUCCESS;

  if (RadiativeTransfer < 1) 
    return SUCCESS;

  if (RadiativeTransfer > 0 && GridRank < 3) {
    fprintf(stderr, "Grid_TransportPhotonPackage: failed\n");
    fprintf(stderr, "Grid_TransportPhotonPackage: "
	    "Transfer in less than 3D is not implemented.\n");
    ENZO_FAIL("");
  }

  if (PhotonPackages->NextPackage == NULL)
    return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stdout, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  int RPresNum1, RPresNum2, RPresNum3;
  if (RadiationPressure)
    IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3);

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    fprintf(stdout, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  if (DEBUG) fprintf(stdout,"TransportPhotonPackage: initialize fields.\n");
  if (DEBUG) fprintf(stdout,"TransportPhotonPackage: %"ISYM" %"ISYM" .\n",
		     GridStartIndex[0], GridEndIndex[0]);

  int i, count, dcount, tcount, pcount, trcount;
  PhotonPackageEntry *PP, *FPP, *SavedPP, *PausedPP;
  PP = PhotonPackages;

  if (DEBUG) {
    count = 0;
    while ((PP->NextPackage) != NULL) { 
      count++;
      PP=PP->NextPackage;
    }
    fprintf(stdout, "TransportPhotonPackage: done initializing.\n");
    fprintf(stdout, "G%d: counted %"ISYM" packages\n", GridNum, count);
  }

  /* If requested, make vertex centered field (only when it doesn't
     exist ... see inside routine). */

  if (RadiativeTransferInterpolateField)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (FieldsToInterpolate[i] == TRUE)
	if (this->ComputeVertexCenteredField(i) == FAIL) {
	  fprintf(stderr, "Error in grid->ComputeVertexCenteredField "
		  "(field %"ISYM").\n", i);
	  ENZO_FAIL("");
	}

  PP = PhotonPackages;

  while (PP != NULL) {

  pcount = dcount = tcount = trcount = 0;

  ListOfPhotonsToMove *LocalMoveList;
  int NumberOfThreads = NumberOfCores / NumberOfProcessors;
#pragma omp parallel private(LocalMoveList,SavedPP,i)	\
  reduction(+:pcount,dcount,tcount,trcount)
  {

    int tpcount, nph, start_num, end_num, thread_num;
    PhotonPackageEntry *ThreadPP0 = new PhotonPackageEntry;
    PhotonPackageEntry *TempPP, *ThreadPP;
    LocalMoveList = new ListOfPhotonsToMove;
    LocalMoveList->NextPackageToMove = NULL;

#ifdef _OPENMP
    thread_num = omp_get_thread_num();
#else
    thread_num = 0;
#endif /* _OPENMP */

    /* Manually split up photon linked list to divide the work among
       cores */

    // Count number of photons in list
    nph = 0;
    TempPP = PhotonPackages;
    while (TempPP->NextPackage != NULL) {
      TempPP = TempPP->NextPackage;
      nph++;
    }

    if (NumberOfThreads > 1 && nph >= NumberOfThreads) {

      start_num = thread_num * nph / NumberOfThreads;
      if (thread_num < NumberOfThreads-1)
	end_num = (thread_num+1) * nph / NumberOfThreads;
      else
	end_num = nph;

      // Find the head of the list on this thread (i=start_num)
      ThreadPP0->PreviousPackage = NULL;
      TempPP = PhotonPackages->NextPackage;
      for (i = 0; i < start_num; i++)
	TempPP = TempPP->NextPackage;

      /* Find the tail of the tail and mark the NextPackage as NULL,
	 so the thread only does its share of work.  We'll piece
	 together the thread lists after the computation. */

      ThreadPP = TempPP;
      ThreadPP->PreviousPackage = ThreadPP0;
      ThreadPP0->NextPackage = ThreadPP;

      for (i = start_num; i < end_num-1; i++)
	TempPP = TempPP->NextPackage;

      /* BARRIER: We don't want to mark the tail node's NextPackage to
	 NULL until every thread has found their head node.
	 Otherwise, it'll run into the NULL NextPackage pointer, set
	 from another thread. */

#pragma omp barrier
      TempPP->NextPackage = NULL;

      // Detach head node from list
#pragma omp single
      PhotonPackages->NextPackage = NULL;
      
    } else {
      ThreadPP = PhotonPackages->NextPackage;
    }

    if (DEBUG>1) {
      TempPP = ThreadPP;
      tpcount = 0;
      printf("T%d: ThreadPP0 = %x\n", thread_num, ThreadPP0);
      while (TempPP != NULL) {
	printf("T%d: photon %d/%d: %x %x %x\n", thread_num, 
	       tpcount, nph, TempPP->PreviousPackage,
	       TempPP, TempPP->NextPackage);
	tpcount++;
	TempPP = TempPP->NextPackage;
      }
    }
    
    count = 0;
    //PP = PhotonPackages->NextPackage;
    FPP = this->FinishedPhotonPackages;
    PausedPP = this->PausedPhotonPackages;
  
    int AdvancePhotonPointer;
    int DeleteMe, DeltaLevel, PauseMe;
    grid *MoveToGrid;

    const float clight = 2.9979e10;
    float LightCrossingTime = 1.7320508 * (LengthUnits/TimeUnits) /
      (clight * RadiativeTransferPropagationSpeedFraction);  // sqrt(3)=1.73
    FLOAT EndTime;
    if (RadiativeTransferAdaptiveTimestep)
      EndTime = PhotonTime+LightCrossingTime;
    else
      EndTime = PhotonTime+dtPhoton-ROUNDOFF;

    while (ThreadPP != NULL) {

      DeleteMe = FALSE;
      PauseMe = FALSE;
      MoveToGrid = NULL;
      AdvancePhotonPointer = TRUE;

      if ((ThreadPP->CurrentTime) < EndTime) {
	WalkPhotonPackage(&ThreadPP,
			  &MoveToGrid, ParentGrid, CurrentGrid, Grids0, nGrids0,
			  DensNum, DeNum, HINum, HeINum, HeIINum, H2INum,
			  kphHINum, gammaNum, kphHeINum, 
			  kphHeIINum, kdissH2INum, RPresNum1,
			  RPresNum2, RPresNum3, DeleteMe, PauseMe, DeltaLevel, 
			  LightCrossingTime,
			  DensityUnits, TemperatureUnits, VelocityUnits, 
			  LengthUnits, TimeUnits);
	tcount++;
      } else {

	/* If all work is finished, store in FinishedPhotonPackages and
	   don't check for work until next timestep */

	SavedPP = PopPhoton(ThreadPP);
	ThreadPP = ThreadPP->NextPackage;
#pragma omp critical
	InsertPhotonAfter(FPP, SavedPP);
	AdvancePhotonPointer = FALSE;

      }

      if (DEBUG > 1) 
	fprintf(stdout, "photon #%"ISYM" %x %x %x\n",
		tcount,  ThreadPP,  PhotonPackages, 
		MoveToGrid); 

      if (PauseMe == TRUE) {
	if (DEBUG > 1) fprintf(stdout, "paused photon %x\n", ThreadPP);
	SavedPP = PopPhoton(ThreadPP);
	//printf("paused photon %x lvl %"ISYM" ipix %"ISYM" CSRC %x leafID %"ISYM"\n",
	//SavedPP, SavedPP->level, SavedPP->ipix, SavedPP->CurrentSource,
	//SavedPP->CurrentSource->LeafID);
	ThreadPP = ThreadPP->NextPackage;
#pragma omp critical
	InsertPhotonAfter(PausedPP, SavedPP);
	AdvancePhotonPointer = FALSE;
	MoveToGrid = NULL;
	pcount++;
      }

      if (DeleteMe == TRUE) {   
	if (DEBUG > 1) fprintf(stdout, "delete photon %x\n", ThreadPP);
	dcount++;
	ThreadPP = DeletePhotonPackage(ThreadPP);
	MoveToGrid = NULL;
      } 

      if (MoveToGrid != NULL) {
	if (DEBUG) {
	  fprintf(stdout, "moving photon from %x to %x\n", 
		  CurrentGrid,  MoveToGrid);
	  fprintf(stdout, "moving photon %x %x %x %x\n", 
		  ThreadPP,  ThreadPP->PreviousPackage, 
		  ThreadPP->NextPackage,  PhotonPackages);
	}
	ListOfPhotonsToMove *NewEntry = new ListOfPhotonsToMove;
	NewEntry->NextPackageToMove = LocalMoveList->NextPackageToMove;
	LocalMoveList->NextPackageToMove = NewEntry;
	NewEntry->PhotonPackage = ThreadPP;
	NewEntry->FromGrid = CurrentGrid;
	NewEntry->ToGrid   = MoveToGrid;
	NewEntry->ToGridNum= MoveToGrid->GetGridID();
	NewEntry->ToLevel  = level + DeltaLevel;
	NewEntry->ToProcessor = MoveToGrid->ReturnProcessorNumber();

	if (NewEntry->ToProcessor >= NumberOfProcessors)
	  printf("TransportPH(P%"ISYM" :: G%"ISYM"): WARNING BAD TO_PROC -- P%"ISYM"->P%"ISYM".\n",
		 MyProcessorNumber, GridNum, ProcessorNumber, 
		 NewEntry->ToProcessor);

	if (ThreadPP->PreviousPackage != NULL) 
	  ThreadPP->PreviousPackage->NextPackage = ThreadPP->NextPackage;
	if (ThreadPP->NextPackage != NULL) 
	  ThreadPP->NextPackage->PreviousPackage = ThreadPP->PreviousPackage;
	trcount++;
      } // ENDIF MoveToGrid

      if (AdvancePhotonPointer == TRUE)
	ThreadPP = ThreadPP->NextPackage;

    } // ENDWHILE photons

    /* After all of the photons are traced, combine all of the list of
       photons to move from each thread into one list */

#pragma omp critical
    {
      ListOfPhotonsToMove **TempL = PhotonsToMove;
      // Find the end of the list and link it there.
      while ((*TempL)->NextPackageToMove != NULL)
	(*TempL) = (*TempL)->NextPackageToMove;
      (*TempL)->NextPackageToMove = LocalMoveList->NextPackageToMove;
    } // END omp critical

    delete LocalMoveList;
    delete ThreadPP0;

  } // END omp parallel

  //  if (DEBUG)
    fprintf(stdout, "grid::TransportPhotonPackage: "
	    "transported %"ISYM" deleted %"ISYM" paused %"ISYM"\n",
	    tcount, dcount, pcount);

  // After threaded calls, reset pointers to the head
  PP = PhotonPackages->NextPackage;
  PausedPP = PausedPhotonPackages;

  // Merge "paused" photons only when all photons have been transported
  if (PP == NULL && PausedPP->NextPackage != NULL) {
    if (this->MergePausedPhotonPackages() == FAIL) {
      fprintf(stderr, "Error in grid::MergePausedPhotonPackages.\n");
      ENZO_FAIL("");
    }
    // Reset temp pointers
    PP = PhotonPackages->NextPackage;
    //FPP = this->FinishedPhotonPackages;
    PausedPP = this->PausedPhotonPackages;
    this->PausedPhotonPackages->NextPackage = NULL;
    this->PausedPhotonPackages->PreviousPackage = NULL;
  } // ENDIF merge photons

  } // ENDWHILE PP != NULL

  NumberOfPhotonPackages -= dcount;

  /* For safety, clean up paused photon list */

#ifdef UNUSED
  if (PausedPhotonPackages->NextPackage != NULL) {
    PausedPP = PausedPhotonPackages->NextPackage;
    while (PausedPP != NULL) {
      PausedPP = DeletePhotonPackage(PausedPP);
      PausedPP = PausedPP->NextPackage;
    }
    PausedPhotonPackages->NextPackage = NULL;
    PausedPhotonPackages->PreviousPackage = NULL;
  }
#endif

  return SUCCESS;
}
