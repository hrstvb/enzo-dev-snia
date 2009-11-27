/***********************************************************************
/
/  OUTPUT ESCAPE FRACTIONS TO FILE
/
/  written by: John Wise
/  date:       November 2007
/  modified1:
/                
/
/  PURPOSE:
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

int OutputEscapeFraction(void)
{

  /* If we're keeping track of photon escape fractions on multiple
     processors, collect photon counts from all processors */

  int i;
  FILE *fptr;

  if (RadiativeTransferPhotonEscapeRadius > 0) {
    CommunicationSumValues(EscapedPhotonCount, 4);
    if (MyProcessorNumber == ROOT_PROCESSOR) {

      /* Open f_esc file for writing */

      if (TotalEscapedPhotonCount[0] <= 0) {
	if ((fptr = fopen(PhotonEscapeFilename, "w")) == NULL) {
	  fprintf(stderr, "Error opening file %s\n", PhotonEscapeFilename);
	  ENZO_FAIL("");
	}
	fprintf(fptr, 
		"# Time TotalPhotons fesc(0.5rvir) fesc(rvir) fesc(2rvir)\n");
      } else {
	if ((fptr = fopen(PhotonEscapeFilename, "a")) == NULL) {
	  fprintf(stderr, "Error opening file %s\n", PhotonEscapeFilename);
	  ENZO_FAIL("");
	}
      }

      for (i = 0; i < 4; i++)
	TotalEscapedPhotonCount[i] += EscapedPhotonCount[i];

      fprintf(fptr, "%"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
	      PhotonTime, TotalEscapedPhotonCount[0],
	      TotalEscapedPhotonCount[1] / TotalEscapedPhotonCount[0], 
	      TotalEscapedPhotonCount[2] / TotalEscapedPhotonCount[0], 
	      TotalEscapedPhotonCount[3] / TotalEscapedPhotonCount[0]);

      fclose(fptr);

    } // ENDIF ROOT_PROCESSOR

  } // ENDIF RTPhotonEscapeRadius

  return SUCCESS;

}
