/***********
 * ComputeSphericalGravityPotential
 * dcollins.  October 16 2018.  14:58.
 *
 * *********/
#include "preincludes.h"
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "EnzoTiming.h"
#include "performance.h"
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

int CommunicationBroadcastValues(FLOAT *Values, int Number, int BroadcastProcessor);

int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[]){
    if ( SphericalGravity == 0 )
        return SUCCESS;
    fprintf(stderr,"CLOWN calling spherical gravity.  Neat \n");
    if ( SphericalGravityBinNumber < 0 && SphericalGravityBinSize < 0 ){
        fprintf(stderr,"FAILURE  SphericalGravityBinNumber = %"ISYM" && SphericalGravityBinSize = %"FSYM"\n",
         SphericalGravityBinNumber, SphericalGravityBinSize );
        ENZO_FAIL("Improper Bin Size and Number, SphericalGravityComputePotential");
    }
    if ( SphericalGravityBinSize < 0 ){
        SphericalGravityBinSize = ( SphericalGravityOuterRadius - SphericalGravityInnerRadius )/SphericalGravityBinNumber;
    }
    if ( SphericalGravityBinNumber < 0 ){
        SphericalGravityBinNumber = int( SphericalGravityOuterRadius - SphericalGravityInnerRadius )/SphericalGravityBinSize;
    }

    fprintf(stderr,"CLOWN: SphericalGravity Nbins %d of size %0.2e\n",SphericalGravityBinNumber, SphericalGravityBinSize);
    if ( SphericalGravityMassInterior != NULL ){
        delete [] SphericalGravityMassInterior;
    }
    if ( SphericalGravityMassShell != NULL ){
        delete [] SphericalGravityMassShell;
    }
    SphericalGravityMassInterior = new float[SphericalGravityBinNumber];
    SphericalGravityMassShell = new float[SphericalGravityBinNumber];

    LevelHierarchyEntry *Temp;
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->SphericalGravityAddMassToShell();
      Temp = Temp->NextGridThisLevel;
    }
    for(int i=0;i<SphericalGravityBinNumber;i++){
        fprintf(stderr,"CLOWN shell pre NP%d %d %0.2e\n",MyProcessorNumber,i, SphericalGravityMassShell[i]);
    }


    //It might be better to use MPI_Alltoallv; see CommunicationShareParticles.C
    //That takes some more setup, so in the words of Mike, "make it work then make it work fast."
    CommunicationSumValues(SphericalGravityMassShell,SphericalGravityBinNumber);
    CommunicationBroadcastValues(SphericalGravityMassShell,SphericalGravityBinNumber,ROOT_PROCESSOR);

    for(int i=0;i<SphericalGravityBinNumber;i++){
        fprintf(stderr,"CLOWN shell post NP%d %d %0.2e\n",MyProcessorNumber,i, SphericalGravityMassShell[i]);
    }

    return SUCCESS;
}
