/* 
 * SphericalGravityAddMassToShell
 * dcollins.  Oct. 16 2018. 15:15.
 */
#include <stdio.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);
int grid::SphericalGravityAddMassToShell(){
    if ( ProcessorNumber != MyProcessorNumber )
        return SUCCESS;
    fprintf(stderr,"CLOWN in the grid whatzit\n");
    for ( int i=0; i<SphericalGravityBinNumber;i++){
        SphericalGravityMassShell[i] = 0.1;
    }
    return SUCCESS;
}
