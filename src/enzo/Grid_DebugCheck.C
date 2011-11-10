/***********************************************************************
/
/  GRID CLASS (CHECKS FOR NANS)
/
/  written by: Greg Bryan
/  date:       1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
 
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
// Turn the first on to check the gas and dark matter data for nan's.
// Turn the second on just to print the message (i.e. which routine
//     has called it.
 
#define DEBUG_CHECK_OFF
#define TRACE_OFF
 
 
 
 
int grid::DebugCheck(char *message)
{
 
#ifdef TRACE_ON
 
  if (ProcessorNumber == MyProcessorNumber)
    fprintf(stderr, "P(%"ISYM"): %s\n", MyProcessorNumber, message);
 
#endif /* TRACE_ON */
 
#ifdef DEBUG_CHECK_ON
 
  // Return if this doesn't concern us
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  // Set this to zero so the compiler doesn't optimize everything away
 
  int ThisIsZero = GridStartIndex[0] - DEFAULT_GHOST_ZONES, size = 1,
      dim, k1, k2;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  for (k1 = 0; k1 < NumberOfBaryonFields; k1++)
    for (k2 = 0; k2 < size; k2++)
      if (BaryonField[k1][k2+ThisIsZero] != BaryonField[k1][k2]) {
	fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): gas %"ISYM" %"ISYM" %"GSYM"\n", message,
		MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	exit(EXIT_FAILURE);
      }
 
  for (k1 = 0; k1 < GridRank; k1++)
    for (k2 = 0; k2 < NumberOfParticles; k2++)
      if (ParticlePosition[k1][k2] != ParticlePosition[k1][k2+ThisIsZero] ||
	  ParticleVelocity[k1][k2] != ParticleVelocity[k1][k2+ThisIsZero]  ) {
	fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): dm %"ISYM" (%"ISYM"/%"ISYM") %"GSYM" %"GSYM"\n",
		message, MyProcessorNumber, k1, k2, NumberOfParticles,
		ParticlePosition[k1][k2], ParticleVelocity[k1][k2]);
	exit(EXIT_FAILURE);
      }
#ifdef MHDCT 
  int HasProblem=FALSE, FailHard=FALSE;
  if( useMHDCT ){

    
    for(field=0;field<3;field++){
      
      
      for(k=0;k<MagneticDims[field][2];k++)
	for(j=0;j<MagneticDims[field][1];j++)
	  for(i=0;i<MagneticDims[field][0];i++){

      if(field == 0)
	return i+MagneticDims[0][0]*(j+MagneticDims[0][1]*k);

      if(field == 1)
	return i+MagneticDims[1][0]*(j+MagneticDims[1][1]*k);
      if(field == 2)
	  return i+MagneticDims[2][0]*(j+MagneticDims[2][1]*k);
	    
	    if( HasProblem == FALSE ){
	      if( isnan(MagneticField[field][index]) ){
		fprintf(stderr, "DebugCheck[%s](Proc %d nan): MagneticField[%d][%d,%d,%d]= %g Left = %13g %13g %13g Width %13g\n\n", message,
			MyProcessorNumber, field, i,j,k, MagneticField[field][index],
			GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2], CellWidth[0][0]);
		HasProblem=TRUE;
	      }
	      if( isinf(MagneticField[field][index]) ){
		fprintf(stderr, "\nDebugCheck[%s](Proc %d inf): MagneticField[%d][%d,%d,%d]= %g Left = %13g %13g %13g Width %13g\n", message,
			MyProcessorNumber, field, i,j,k, MagneticField[field][index],
			GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2], CellWidth[0][0]);
		HasProblem=TRUE;
	      }
	    }
	  }//i,j,k
      
      if( ElectricField[field] != NULL )
	for(k=0;k<ElectricDims[field][2];k++)
	  for(j=0;j<ElectricDims[field][1];j++)
	    for(i=0;i<ElectricDims[field][0];i++){
	      index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
	  if( (FailHard == TRUE && HasProblem == FALSE ) || FailHard == FALSE ) {
	      if (ElectricField[field][index+ThisIsZero] != ElectricField[field][index]) {
		fprintf(stderr, "DebugCheck[%s](Proc %d): ElectricField[%d][%d,%d,%d]= %g\n", message, 
			MyProcessorNumber, field, i,j,k, ElectricField[field][index]);
		HasProblem=TRUE;
	      }
	  }
	    }//i,j,k
      
    }//field

  }//MHD

#endif //MHDCT
#if 0		
  if (NumberOfBaryonFields > 0)
    for (k1 = 0; k1 < 2+DualEnergyFormalism; k1++)
      for (k2 = 0; k2 < size; k2++)
	if (BaryonField[k1][k2+ThisIsZero] <= 0) {
	  fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): <0 %"ISYM" %"ISYM" %"GSYM"\n", message,
		  MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	  exit(EXIT_FAILURE);
	}
#endif
 
#endif /* DEBUG_CHECK_ON */
		
  return SUCCESS;
}
