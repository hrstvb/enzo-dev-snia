
//
// ExtraOutput
// 
// Written by: David Collins 
// date      : May 25, 2011.  3:36 pm.  Cold in my office.
// 
// Purpose   : To provide a horrific amount of data.  Calls to ExtraOutput(int) are sprinkled throughout the code,
//             which helps debug events throughout the evolution.
#include "preincludes.h"
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "performance.h"
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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1, int CheckpointDump = FALSE);
static int output_number[] = {0,0,0,0,0,0,0,0,0,0};
int ExtraOutput(int output_flag, LevelHierarchyEntry *LevelArray[],TopGridData *MetaData, int level, ExternalBoundary *Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
        ){
    int WriteOut = FALSE;
    for( int i=0; i<MAX_EXTRA_OUTPUTS; i++){
        fprintf(stderr,"CRACKER lies ExtraOutputs[%d]=%d, output_flag = %d, %d\n",i,ExtraOutputs[i],output_flag, output_flag - ExtraOutputs[i]);
        if( output_flag == ExtraOutputs[i]){
            WriteOut = TRUE;
            fprintf(stderr,"FUCK THIS SHIT\n");
        }
    }
        fprintf(stderr,"DERPx %s %d\n",MetaData->ExtraDumpName,TRUE);
    if( WriteOut ){
        LevelHierarchyEntry *Temp2 = LevelArray[0];
        while (Temp2->NextGridThisLevel != NULL)
          Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
        //#ifdef USE_HDF5_GROUPS
        sprintf(MetaData->ExtraDumpName,"ExtraDump%02d",output_flag);
        sprintf(MetaData->ExtraDumpDir,"ED%02d",output_flag);
        fprintf(stderr,"DERP2 %s\n",MetaData->ExtraDumpName);
        if (Group_WriteAllData(MetaData->ExtraDumpName, output_number[output_flag]++,
                   Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef TRANSFER
                   ImplicitSolver,
#endif
                   LevelArray[level]->GridData->ReturnTime(), FALSE) == FAIL) {
                ENZO_FAIL("Error in Group_WriteAllData.");
        }
    }
  return SUCCESS;
}
