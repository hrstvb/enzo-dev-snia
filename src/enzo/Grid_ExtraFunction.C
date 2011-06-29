//
// ExtraFunction
// Generic grid member function for debugging.
//

#include <stdio.h>
#include "ErrorExceptions.h"
#include "math.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
static float LastB[] = {0.0,0.0,0.0};
int grid::ExtraFunction(char * message){
    if( MyProcessorNumber != ProcessorNumber ){
        return SUCCESS;
    }
    int size=0;
    float total=0;
    fprintf(stderr,"MagField %s",message);
    for(int field=0;field<3;field++){
        size =MagneticDims[field][0]*MagneticDims[field][1]*MagneticDims[field][2];
        total=1;

        for(int index=0;index<size; index++){
            total += fabs(MagneticField[field][index]);
        }
        fprintf(stderr," %"FSYM,total-LastB[field]);
        LastB[field]=total;
    }
    fprintf(stderr,"\n");

    return SUCCESS;

}
