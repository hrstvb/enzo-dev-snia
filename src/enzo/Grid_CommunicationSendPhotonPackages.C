#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (SEND PHOTONS FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE: 
/
/    NOTE: We assume all the from grids are at the same level!
/    NOTE: Modelled after GBs Grid_CommunicationSendParticles.C
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "GroupPhotonList.h"
#include "Parallel.h"
#include "CommunicationUtilities.h"

using namespace Parallel;

void my_exit(int status);
int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);

int grid::CommunicationSendPhotonPackages(grid *ToGrid, int ToProcessor,
					  int ToNumber, int FromNumber, 
					  PhotonPackageEntry **ToPP,
					  int CommunicationIndex)
{

  int index, dim, temp_int;
  PhotonPackageEntry *PP;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (FromNumber == 0) 
    return SUCCESS;

  /* Allocate MPI buffer */

  const int CommType = 15;
  MPIBuffer *mbuffer = NULL;
  if (Parallel::CommunicationDirection != COMMUNICATION_RECEIVE) {
    int iarg[] = {ToNumber, FromNumber, 0};
    mbuffer = new MPIBuffer(this, ToGrid, CommType, MPI_PHOTON_TAG,
			    NULL, NULL, NULL, iarg);
  } else {
    mbuffer = GetMPIBuffer(CommunicationIndex);  // Grab from list.
  }
 
  /* Allocate memory */

  PhotonBuffer *buffer = NULL;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (PhotonBuffer *) mbuffer->ReturnBuffer();
  else
#endif /* USE_MPI */
    buffer = new PhotonBuffer[FromNumber];

  if (DEBUG)
    printf("SendPhotonPackages(%"ISYM"): Sending %"ISYM" photon packages "
	   "from P%"ISYM"->P%"ISYM".\n",
	   MyProcessorNumber, FromNumber, ProcessorNumber, ToProcessor);

  /* If this is from processor, pack photons */

  if (MyProcessorNumber == ProcessorNumber) {
    index = 0;
    PP = PhotonPackages->NextPackage;

    while (PP != NULL) {
      buffer[index].Photons		 = PP->Photons;
      buffer[index].Type		 = PP->Type;
      buffer[index].Energy		 = PP->Energy;
      buffer[index].EmissionTimeInterval = PP->EmissionTimeInterval;
      buffer[index].EmissionTime	 = PP->EmissionTime;
      buffer[index].CurrentTime          = PP->CurrentTime;
      buffer[index].ColumnDensity        = PP->ColumnDensity;
      buffer[index].CrossSection         = PP->CrossSection;
      buffer[index].Radius		 = PP->Radius;
      buffer[index].ipix		 = PP->ipix;
      buffer[index].level		 = PP->level;
      for (dim = 0; dim < GridRank; dim++)
	buffer[index].SourcePosition[dim] = PP->SourcePosition[dim];
      buffer[index].SourcePositionDiff   = PP->SourcePositionDiff;

      if (PP->CurrentSource != NULL)
	buffer[index].SuperSourceID = PP->CurrentSource->LeafID;
      else
	buffer[index].SuperSourceID = -1;
      
      if (PP->CurrentTime < 0 || PP->CurrentTime > 1e10) {
	ENZO_VFAIL("CTPhotons[0][P%"ISYM"->P%"ISYM"]: "
		"(%"ISYM" of %"ISYM") Bad photon time %"GSYM"\n",
		ProcessorNumber, ToProcessor, index, NumberOfPhotonPackages, 
		PP->CurrentTime)
      }

      // Next photon
      PP = PP->NextPackage;
      index++;

    }  /* ENDWHILE PP != NULL */

    if (DEBUG)
      printf("CommSendPhotons(P%"ISYM"): Counted %"ISYM" photons.\n", MyProcessorNumber,
	     index);

    /* Now that we're done packing the photons, delete them */
    
    PP = PhotonPackages->NextPackage;
    while (PP != NULL) {
      PP = DeletePhotonPackage(PP);
      PP = PP->NextPackage;
    }

    /* Check if we packed all of the photons */

    if (index != FromNumber) {
      fprintf(stdout, "CommSendPhotons WARNING: Counted %"ISYM" photon packages, but"
	      " FromNumber = %"ISYM"\n", index, FromNumber);
      FromNumber = min(index, FromNumber);
      fprintf(stdout, "CommSendPhotons: Correcting FromNumber to %"ISYM"\n", 
	      FromNumber);
      //ENZO_FAIL("Photon package mismatch!\n");
    }

  } /* ENDIF PackPhotons */

#ifdef USE_MPI

  /* Send buffer if the processor numbers aren't identical */

  if (ProcessorNumber != ToProcessor) {

    if (CommunicationDirection != COMMUNICATION_RECEIVE)
      mbuffer->FillBuffer(MPI_PhotonBuffer, FromNumber, buffer);

    if (MyProcessorNumber == ProcessorNumber) {
      if (DEBUG)
	printf("PhotonSend(P%"ISYM"): Sending %"ISYM" photons to processor %"ISYM".\n",
	       MyProcessorNumber, FromNumber, ToProcessor);
      mbuffer->SendBuffer(ToProcessor);
    }

    if (MyProcessorNumber == ToProcessor) {

      if (DEBUG) 
	printf("PhotonSend(P%"ISYM"): Receiving %"ISYM" photons from processor %"ISYM".\n",
	       MyProcessorNumber, FromNumber, ProcessorNumber);

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	mbuffer->IRecvBuffer(ProcessorNumber);
      } // ENDIF post receive

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	mbuffer->RecvBuffer(ProcessorNumber);

    } /* ENDIF (MyProcessorNumber == ToProcessor) */

  } /* ENDIF (ProcessorNumber != ToProcessor) */
#endif /* USE_MPI */

  /* If this is the to processor, unpack fields */

  PhotonPackageEntry *NewPP;

  if (MyProcessorNumber == ToProcessor && 
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    for (index = 0; index < FromNumber; index++) {

      NewPP = new PhotonPackageEntry;

      // insert pointer in list after ToPP
      NewPP->NextPackage = (*ToPP)->NextPackage;
      (*ToPP)->NextPackage = NewPP;
      NewPP->PreviousPackage = *ToPP;
      if (NewPP->NextPackage != NULL)
	NewPP->NextPackage->PreviousPackage = NewPP;

      // Unpack photons
      NewPP->Photons		  = buffer[index].Photons;
      NewPP->Type		  = buffer[index].Type;
      NewPP->EmissionTimeInterval = buffer[index].EmissionTimeInterval;
      NewPP->EmissionTime	  = buffer[index].EmissionTime;
      NewPP->CurrentTime	  = buffer[index].CurrentTime;
      NewPP->ColumnDensity	  = buffer[index].ColumnDensity;
      NewPP->CrossSection	  = buffer[index].CrossSection;
      NewPP->Radius		  = buffer[index].Radius;
      NewPP->ipix		  = buffer[index].ipix;
      NewPP->level		  = buffer[index].level;
      NewPP->Energy		  = buffer[index].Energy;
      for (int dim = 0; dim < GridRank; dim++) 
	NewPP->SourcePosition[dim]  = buffer[index].SourcePosition[dim];
      NewPP->SourcePositionDiff   = buffer[index].SourcePositionDiff;

      if (NewPP->CurrentTime < 0 || NewPP->CurrentTime > 1e10) {
	ENZO_VFAIL("CTPhotons[1][P%"ISYM"->P%"ISYM"]: "
		"(%"ISYM" of %"ISYM") Bad photon time %"GSYM"\n",
		ProcessorNumber, ToProcessor, index, FromNumber, 
		NewPP->CurrentTime)
      }

      if (RadiativeTransferSourceClustering) {
	if (FindSuperSource(&NewPP, buffer[index].SuperSourceID) == FAIL) {
	  ENZO_FAIL("Error in FindSuperSource.\n");

	}
      } else
	NewPP->CurrentSource = NULL;

      // Move pointer to the next photon
      //      *ToPP = (*ToPP)->NextPackage;

    } /* ENDFOR index */

    /* Only delete the buffer if we're in receive mode (in send mode
       it will be deleted by CommunicationBufferedSend and if we're in
       post-receive mode then it will be deleted when we get to
       receive-mode). */

    delete [] buffer;

  }  /* ENDIF (MyProcessorNumber == ToProcessor) */

  return SUCCESS;

}
