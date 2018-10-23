/* 
 * SphericalGravityAddMassToShell
 * dcollins.  Oct. 16 2018. 15:15.
 */
#include <stdio.h>
#include "hdf5.h"
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

    FLOAT x,y,z,r;
    FLOAT CellVolume = CellWidth[0][0];
    if( GridRank > 1 ) CellVolume *= CellWidth[1][0];
    if( GridRank > 2 ) CellVolume *= CellWidth[2][0];

    int rbin, index;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
                                   TENum, B1Num, B2Num, B3Num);

    for( int k=GridStartIndex[2]; k<GridEndIndex[2];k++)
    for( int j=GridStartIndex[1]; j<GridEndIndex[1];j++)
    for( int i=GridStartIndex[0]; i<GridEndIndex[0];i++){
        index = i + GridDimension[0]*( j + GridDimension[1]*k);
        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphericalGravityCenter[0];
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SphericalGravityCenter[1];
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SphericalGravityCenter[2];
        r=sqrt( x*x+y*y+z*z );
        if ( r >= SphericalGravityInnerRadius && r <= SphericalGravityOuterRadius ){
            rbin = int( (r-SphericalGravityInnerRadius) / SphericalGravityBinSize );
            SphericalGravityMassShell[rbin] += BaryonField[DensNum][index] * CellVolume;
            SphericalGravityBinCount[rbin] += 1;
        }
    }
    return SUCCESS;
}

