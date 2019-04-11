/***********************************************************************
 /
 /  GRID CLASS (ALLOCATE SPACE FOR GRIDS -- DIMS, ETC MUST BE SET)
 /
 /  written by: Greg Bryan
 /  date:       July, 1995
 /  modified1:
 /
 /  PURPOSE:
 /
 ************************************************************************/

//
//  Allocate room for the grids.
//
#include <stdio.h>
#include <stdlib.h>

#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::AllocateGrids()
{

	/* Compute grid size. */

//  int size = 1,i,field;
//  for (int dim = 0; dim < GridRank; dim++)
//    size *= GridDimension[dim];
//
//  /* Allocate room and clear it. */
//
//  for (field = 0; field < NumberOfBaryonFields; field++) {
//    BaryonField[field]    = new float[size];
//    for (i = 0; i < size; i++)
//      BaryonField[field][i] = 0.0;
//  }
//
//  if(UseMHDCT){
//    for(field=0;field<3;field++){
//
//
//      MagneticField[field] = new float[MagneticSize[field]];
//      ElectricField[field] = new float[ElectricSize[field]];
//      AvgElectricField[field] = new float[ ElectricSize[field] ];
//
//      for(i=0; i<ElectricSize[field]; i++) ElectricField[field][i] = 0.0;
//      for(i=0; i<MagneticSize[field]; i++) MagneticField[field][i] = 0.0;
//
//      MHDParentTemp[field] = NULL;
//    }
//
//  } // if(UseMHDCT)
	size_t size = GetGridSize();
	for(int field = 0; field < NumberOfBaryonFields; field++)
		arr_newset(BaryonField + field, size, 0);

	if(UseMHDCT)
	{
		for(int dim = 0; dim < 3; dim++)
		{
			arr_newset(MagneticField + dim, MagneticSize[dim], 0);
			arr_newset(ElectricField + dim, ElectricSize[dim], 0);
			arr_newset(AvgElectricField + dim, ElectricSize[dim], 0);
			MHDParentTemp[dim] = NULL;
		}
	}

	if(RandomForcing == TRUE)
		for(int dim = 0; dim < GridRank; dim++)
			arr_newset(RandomForcingField + dim, size, 1.0);
}
