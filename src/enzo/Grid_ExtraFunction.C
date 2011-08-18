//
// ExtraFunction
// Generic grid member function for debugging.
// 
// A stub function that accepts only a string as an argument.
//
// One possible use: checking for mass conservation errors in the code.  By filling
// the body with something that loops over zones and totals mass, one can then look for leaks
// by calling this code from various points.   
//
// All code in the body of the function should be considered temporary,
// even if checked into the repo.
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
  return SUCCESS;
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
