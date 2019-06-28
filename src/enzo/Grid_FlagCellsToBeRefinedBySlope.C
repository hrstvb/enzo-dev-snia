/***********************************************************************
 /
 /  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE)
 /
 /  written by: Greg Bryan
 /  date:       November, 1994
 /  modified1:  Alexei Kritsuk, Feb. 2005. Fixed a bug by adding the second
 /                              fabs() at around L107. Added a switch on
 /                              ProblemType.
 /  modified2:  Brian O'Shea, May 2006.  Changed problemtype 30 to just refine
 /                              on Density and TotalEnergy (as suggested by
 /                              Alexei Kritsuk).  This is for supernovae.
 /
 /  PURPOSE:
 /
 /  RETURNS:
 /    number of flagged cells, or -1 on failure
 /
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "DebugMacros.h"
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::FlagCellsToBeRefinedBySlope()
{
	/* Return if this grid is not on this processor. */
	if(MyProcessorNumber != ProcessorNumber)
		return SUCCESS;

	/* error check */
	if(FlaggingField == NULL)
	{
		fprintf(stderr, "Flagging Field is undefined.\n");
		return -1;
	}

	float MinimumSlopeForRefinementThis = MinimumSlopeForRefinement[0]; // default
	bool doAll = SlopeFlaggingFields[0] == INT_UNDEFINED;
	int Start[MAX_DIMENSION], End[MAX_DIMENSION];
	arr_set(Start, MAX_DIMENSION, 0);
	arr_set(End, MAX_DIMENSION, 0);
	arr_cpy(Start, GridStartIndex, GridRank);
	arr_cpy(End, GridEndIndex, GridRank);

	//////////////////////////////////////////////////////////////////////
	//
	// BARYON FIELDS
	//
	/* Some problems do not need to check all the fields, in particular
	 the velocity components */
	int indexInFlaggingList, thisFieldType;
	bool isBaryonField;
	int nBaryonFields = NumberOfBaryonFields;

	// Override NumberOfFields for some specific problems:
	if(ProblemType == 6)
		nBaryonFields = 1; // Implosion (AK)
	if(ProblemType == 7)
		nBaryonFields = 2; // SedovBlast (AK)
	if(ProblemType == 11)
		nBaryonFields = 2; // RadiatingShock (BWO)
	if(ProblemType == 30)
		nBaryonFields = 2; // Cosmology (BWO 23 May 2006)

	//	for(dim = 0; dim < GridRank; dim++) {
	//		if(GridDimension[dim] <= 1)
	//			continue;
	//
	//		for(int field = 0; field < nBaryonFields; field++) {
	//			bool doField = false;
	//			float MinimumSlopeForRefinementThis;
	//			if(SlopeFlaggingFields[0] == INT_UNDEFINED) {
	//				MinimumSlopeForRefinementThis = MinimumSlopeForRefinement[0];
	//				doField = true;
	//			}
	//			else {
	//				for(int g = 0; g < MAX_FLAGGING_METHODS; g++) {
	//					if(SlopeFlaggingFields[g] == FieldType[field]) {
	//						MinimumSlopeForRefinementThis = MinimumSlopeForRefinement[g];
	//						doField = true;
	//					}
	//				}
	//			}
	//
	//			if(doField)
	//			{
	//				...
	//			}

	for(int field = 0; field < nBaryonFields; field++)
	{
		thisFieldType = FieldType[field];
		if(!doAll)
		{
			int indexInFlaggingList = indexof(SlopeFlaggingFields, MAX_FLAGGING_METHODS, thisFieldType);
			if(indexInFlaggingList < 0)
				continue;
			MinimumSlopeForRefinementThis = MinimumSlopeForRefinement[indexInFlaggingList];
		}
		float* fieldData = BaryonField[field];
		bool useRelativeGrowth = thisFieldType != BurnedMolarFraction;

		for(int dim = 0; dim < GridRank; dim++)
		{
			if(GridDimension[dim] <= 1)
				continue;

			long long Offset = getGridStride(dim);
			float* fieldData1 = fieldData - Offset;
			float* fieldData2 = fieldData + Offset;
			for(size_t k = Start[2]; k <= End[2]; k++)
			{
				for(size_t j = Start[1]; j <= End[1]; j++)
				{
					size_t index = ELT(Start[0], j, k);
					for(size_t i = Start[0]; i <= End[0]; index++, i++)
					{
						if(FlaggingField[index])
							continue;

						//*(TempBuffer + index) = 0.5
						//	* fabs((*(BaryonField[field] + index + Offset)
						//			- *(BaryonField[field] + index - Offset))
						//	/ max(fabs(*(BaryonField[field] + index)), tiny_number));
						float u = fabs(*fieldData);
						if(u <= tiny_number)
							FlaggingField[index] = 1;
						else
						{
							float growth = 0.5 * fabs(fieldData2[index] - fieldData1[index]);
							if(useRelativeGrowth)
								growth /= u;
							FlaggingField[index] = growth > MinimumSlopeForRefinementThis;
						}
					} // end loop over i
				} // end loop over j
			} // end loop over k
		} // end loop over dim

		//  for(i = 0; i < GetGridSize(); i++)
		//		FlaggingField[i] += (TempBuffer[i] > MinimumSlopeForRefinementThis) ? 1 : 0;
	} // end loop over fields

	//////////////////////////////////////////////////////////////////////
	//
	// CALCULATED FIELDS
	//
	thisFieldType = BurnedMolarFraction;
	indexInFlaggingList = indexof(SlopeFlaggingFields, MAX_FLAGGING_METHODS, thisFieldType);
	isBaryonField = 0 <= indexof(FieldType, nBaryonFields, thisFieldType);
	if((0 <= indexInFlaggingList || doAll) && !isBaryonField)
	{
		if(0 <= indexInFlaggingList)
			MinimumSlopeForRefinementThis = MinimumSlopeForRefinement[indexInFlaggingList];

		float* rhoData = NULL;
		float* niData = NULL;
		MHD_SNIA_GetFields(&rhoData, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &niData, NULL, NULL);
		if(rhoData != NULL && niData != NULL)
		{
			for(int dim = 0; dim < GridRank; dim++)
			{
				if(GridDimension[dim] <= 1)
					continue;

				long long Offset = getGridStride(dim);
				float* rhoData1 = rhoData - Offset;
				float* rhoData2 = rhoData + Offset;
				float* niData1 = niData - Offset;
				float* niData2 = niData + Offset;
				for(size_t k = Start[2]; k <= End[2]; k++)
				{
					for(size_t j = Start[1]; j <= End[1]; j++)
					{
						size_t index = ELT(Start[0], j, k);
						for(size_t i = Start[0]; i <= End[0]; index++, i++)
						{
							if(FlaggingField[index])
								continue;

							float a = niData1[index]; //rho_ni left
							a = a / (4 * rhoData1[index] - 3 * a); // burned fraction left
							float b = niData2[index]; //rho_ni right
							b = b / (4 * rhoData2[index] - 3 * b); // burned fraction right
							FlaggingField[index] = fabs(b - a) > MinimumSlopeForRefinementThis;
						} // end loop over i
					} // end loop over j
				} // end loop over k
			}  // end loop over dim
		}
	} // end if BurnedMolarFraction

	//Offset *= GridDimension[dim];
	//delete[] TempBuffer;

	/* Count number of flagged Cells. */
	// Note: The count will include all flagged cells whether flagged by this method or before.
	size_t NumberOfFlaggedCells = 0;
	for(size_t i = 0; i < gridSize; i++)
		if(FlaggingField[i] > 0)
			NumberOfFlaggedCells++;

	return NumberOfFlaggedCells;
}

///***********************************************************************
///
///  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE)
///
///  written by: Greg Bryan
///  date:       November, 1994
///  modified1:  Alexei Kritsuk, Feb. 2005. Fixed a bug by adding the second
///                              fabs() at around L107. Added a switch on
///                              ProblemType.
///  modified2:  Brian O'Shea, May 2006.  Changed problemtype 30 to just refine
///                              on Density and TotalEnergy (as suggested by
///                              Alexei Kritsuk).  This is for supernovae.
///
///  PURPOSE:
///
///  RETURNS:
///    number of flagged cells, or -1 on failure
///
//************************************************************************/
//
//#include <stdio.h>
//#include <math.h>
//#include "ErrorExceptions.h"
//#include "macros_and_parameters.h"
//#include "typedefs.h"
//#include "global_data.h"
//#include "Fluxes.h"
//#include "GridList.h"
//#include "ExternalBoundary.h"
//#include "Grid.h"
//
//int grid::FlagCellsToBeRefinedBySlope()
//{
//  /* declarations */
//
//  int i, j, k, index, dim;
//  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
//
//  /* Return if this grid is not on this processor. */
//
//  if (MyProcessorNumber != ProcessorNumber)
//    return SUCCESS;
//
//  /* error check */
//
//  if (FlaggingField == NULL) {
//    fprintf(stderr, "Flagging Field is undefined.\n");
//    return -1;
//  }
//
//  /* Make sure quantities are defined at least to dim 3 */
//
//  for (dim = GridRank; dim < 3; dim++) {
//    GridDimension[dim] = 1;
//    GridStartIndex[dim] = 0;
//    GridEndIndex[dim] = 0;
//  }
//
//  /* loop over all zones */
//
//  for (dim = 0; dim < 3; dim++) {
//    Start[dim] = GridStartIndex[dim];
//    End[dim]   = GridEndIndex[dim];
//  }
//
//  /* compute size */
//
//  int size = 1;
//  for (dim = 0; dim < GridRank; dim++)
//    size *= GridDimension[dim];
//
//  /* allocate a temporary slope field. */
//
//  float *TempBuffer = new float[size];
//
//  /* loop over active dimensions */
//
//  int Offset = 1;
//
//  /* some problems do not need to check all the fields, in particular
//     the velocity components */
//
//  int NumberOfFields = NumberOfBaryonFields;
//
//  // Override NumberOfFields for some specific problems:
//
//  if (ProblemType ==  6) NumberOfFields = 1; // Implosion (AK)
//  if (ProblemType ==  7) NumberOfFields = 2; // SedovBlast (AK)
//  if (ProblemType == 11) NumberOfFields = 2; // RadiatingShock (BWO)
//  if (ProblemType == 30) NumberOfFields = 2; // Cosmology (BWO 23 May 2006)
//
//  for (dim = 0; dim < GridRank; dim++)
//    if (GridDimension[dim] > 1) {
//
//      for (int field = 0; field < NumberOfFields; field++) {
//
//
//	bool doField=false; float MinimumSlopeForRefinementThis;
//	if (SlopeFlaggingFields[0]==INT_UNDEFINED){
//	  MinimumSlopeForRefinementThis=MinimumSlopeForRefinement[0];
//	  doField=true;}
//	else{
//	  for (int g=0; g<MAX_FLAGGING_METHODS; g++){
//	     if (SlopeFlaggingFields[g]==FieldType[field]){
//	       MinimumSlopeForRefinementThis=MinimumSlopeForRefinement[g];
//	       doField=true;
//	     }
//	  }
//	}
//
//
//
//	  if (doField){
//
//
//	/* zero slope */
//
//	for (i = 0; i < size; i++)
//	  *(TempBuffer + i) = 0.0;
//
//	/* compute slope */
//
//	for (k = Start[2]; k <= End[2]; k++)
//	  for (j = Start[1]; j <= End[1]; j++)
//	    for (i = Start[0]; i <= End[0]; i++) {
//	      index = i + j*GridDimension[0] +
//		      k*GridDimension[1]*GridDimension[0];
//	      *(TempBuffer + index) = 0.5*fabs(
//		     (*(BaryonField[field] + index + Offset) -
//		      *(BaryonField[field] + index - Offset)  ) /
//		  max(fabs(*(BaryonField[field] + index) ), tiny_number));
//	    }
//
//	/* flag field based on slope */
//
//	for (i = 0; i < size; i++)
//	  FlaggingField[i] += (TempBuffer[i] > MinimumSlopeForRefinementThis) ? 1 : 0;
//
//	  } // end loop over do condition
//
//      }  // end loop over field
//
//      Offset *= GridDimension[dim];
//
//    }  // end loop over dimension
//
//  /* delete buffer */
//
//  delete [] TempBuffer;
//
//  /* Count number of flagged Cells. */
//
//  int NumberOfFlaggedCells = 0;
//  for (i = 0; i < size; i++)
//    if (FlaggingField[i] > 0)
//      NumberOfFlaggedCells++;
//
//  return NumberOfFlaggedCells;
//
//}
