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
	int debug_level = 0;

	if(debug_level > 0)
		TRACEGF("Allocating %lld baryon fields of size %lld:", NumberOfBaryonFields, size);
	for(int field = 0; field < NumberOfBaryonFields; field++)
	{
		if(debug_level > 1)
		{
			if(DataLabel[field])
				TRACEGF(" ... allocating field type %lld, %s ...", FieldType[field], DataLabel[field]);
			else
				TRACEGF(" ... allocating field type %lld ...", FieldType[field]);
		}
		arr_newset(BaryonField + field, size, 0);
	}

	if(UseMHDCT)
	{
		for(int dim = 0; dim < 3; dim++)
		{
			if(debug_level > 1)
				TRACEGF(" ... allocating MagneticField[%lld] ...", dim);
			arr_newset(MagneticField + dim, MagneticSize[dim], 0);
			if(debug_level > 1)
				TRACEGF(" ... allocating ElectricField[%lld] ...", dim);
			arr_newset(ElectricField + dim, ElectricSize[dim], 0);
			if(debug_level > 1)
				TRACEGF(" ... allocating AvgElectricField[%lld] ...", dim);
			arr_newset(AvgElectricField + dim, ElectricSize[dim], 0);
			MHDParentTemp[dim] = NULL;
		}
	}

	if(RandomForcing == TRUE)
		for(int dim = 0; dim < GridRank; dim++)
			arr_newset(RandomForcingField + dim, size, 1.0);

	if(debug_level > 0)
		TRACEGF("Allocating baryon fields done.");
}
