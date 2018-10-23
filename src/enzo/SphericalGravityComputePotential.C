/***********
 * ComputeSphericalGravityPotential
 * dcollins.  October 16 2018.  14:58.
 * Radially bins mass into shells and interior mass.
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

    if ( SphericalGravityMassInterior != NULL ){
        delete [] SphericalGravityMassInterior;
    }
    if ( SphericalGravityMassShell != NULL ){
        delete [] SphericalGravityMassShell;
    }
    if ( SphericalGravityBinCount != NULL ){
        delete [] SphericalGravityBinCount;
    }
    if ( SphericalGravityBinCenters == NULL ){
        delete [] SphericalGravityBinCenters;
    }
    SphericalGravityMassInterior = new float[SphericalGravityBinNumber];
    SphericalGravityMassShell = new float[SphericalGravityBinNumber];
    SphericalGravityBinCount  = new float[SphericalGravityBinNumber];
    SphericalGravityBinCenters = new FLOAT[SphericalGravityBinNumber];

    for ( int i=0; i<SphericalGravityBinNumber; i++){
        SphericalGravityMassInterior[i] = 0.0;
        SphericalGravityMassShell[i]    = 0.0;
        SphericalGravityBinCount[i]     = 0;
        SphericalGravityBinCenters[i] = SphericalGravityInnerRadius + (i+0.5)*SphericalGravityBinSize;
    }

    LevelHierarchyEntry *Temp;
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->SphericalGravityAddMassToShell();
      Temp = Temp->NextGridThisLevel;
    }


    //It might be better to use MPI_Alltoallv; see CommunicationShareParticles.C
    //That takes some more setup, so in the words of Mike, "make it work then make it work fast."
    CommunicationSumValues(SphericalGravityMassShell,SphericalGravityBinNumber);
    CommunicationBroadcastValues(SphericalGravityMassShell,SphericalGravityBinNumber,ROOT_PROCESSOR);
    //Bin Count is only used for debugging.
    CommunicationSumValues(SphericalGravityBinCount,SphericalGravityBinNumber);

   SphericalGravityMassInterior[0]=SphericalGravityMassShell[0];
   for(int i=1;i<SphericalGravityBinNumber;i++){
       SphericalGravityMassInterior[i]= SphericalGravityMassInterior[i-1]+
                                        SphericalGravityMassShell[i];
   }

    return SUCCESS;
}
int SphericalGravityWritePotential(char * name ) {
  if ( MyProcessorNumber != ROOT_PROCESSOR  || SphericalGravity == 0){
      return SUCCESS;
  }
  if( SphericalGravityMassInterior == NULL ){
      fprintf(stderr,"SphericalGravityWritePotential: No Mass Defined, not writing.\n");
      return SUCCESS;
  }
  
  hid_t       file_id, mass_dset,shell_dset, radius_dset, count_dset, dataspace_id;
  herr_t      status, h5_status, h5_error = -1;
  
  char filename[100];
  sprintf(filename, "%s.SphericalGravity.h5",name);
  
  //file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  if( file_id == -1 ){
      fprintf(stderr,"ERROR IN ERROR: ignore previous warning.  Opening hdf5 file.\n");
      file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
  }
  
  //Create Dataspace 
  hsize_t number[1] = {SphericalGravityBinNumber};
  dataspace_id=H5Screate_simple(1, number, NULL);
  
  //create set
  //                       above, name,      datatype,  shape of data, Something I dont get
  mass_dset = H5Dcreate(file_id, "SphericalGravityMassInterior", HDF5_PREC, dataspace_id, H5P_DEFAULT);
  shell_dset = H5Dcreate(file_id, "SphericalGravityMassShell", HDF5_PREC, dataspace_id, H5P_DEFAULT);
  radius_dset = H5Dcreate(file_id, "SphericalGravityRadius", HDF5_PREC, dataspace_id, H5P_DEFAULT);
  count_dset = H5Dcreate(file_id, "SphericalGravityBinCount", HDF5_PREC, dataspace_id, H5P_DEFAULT);
  
  //Write the Data Set
  //                (set, memory type, mem. space, file space, transfer details, actual data)
   status = H5Dwrite(mass_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityMassInterior);
   status = H5Dwrite(shell_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityMassShell);
   status = H5Dwrite(radius_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinCenters);
   status = H5Dwrite(count_dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, SphericalGravityBinCount);
   
   
   status = H5Sclose(dataspace_id);
   status = H5Dclose(mass_dset);
   status = H5Dclose(shell_dset);
   status = H5Dclose(radius_dset);
   status = H5Dclose(count_dset);
   status = H5Fclose(file_id);
}
