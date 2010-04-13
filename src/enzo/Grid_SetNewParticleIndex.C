/***********************************************************************
/
/  GRID CLASS (SEARCH FOR ALL STAR PARTICLES AND RETURN HOW MANY)
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:  JHK & JHW (2009)
/
/  PURPOSE:
/
/  NOTE:
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
#include "Grid.h"

void grid::SetNewParticleIndex(int &NumberCount1, PINT &NumberCount2)
{
  FILE * fp;
  char fname[MAX_LINE_LENGTH];
  sprintf(fname, "NewStars%03d.out", MyProcessorNumber);
  fp = fopen(fname, "a");
  
  // Figure the grid cell size
  int num = GridDimension[0] - 2*DEFAULT_GHOST_ZONES;
  float fnum = (float)(num);
  float cell = (GridRightEdge[0] - GridLeftEdge[0])/fnum;
  
  int n, abstype;
  for (n = 0; n < NumberOfParticles; n++) 
    if (ParticleNumber[n] == INT_UNDEFINED) {
      abstype = ABS(ParticleType[n]);
      if (abstype == PARTICLE_TYPE_STAR ||
	  (abstype >= PARTICLE_TYPE_SINGLE_STAR &&
	   abstype != PARTICLE_TYPE_MBH))
	ParticleNumber[n] = NumberCount1++ + NumberCount2;
      else 
	ParticleNumber[n] = NumberCount1 + NumberCount2++;

    fprintf(fp, "%d %1.12e %1.12e %1.12e %1.12e %1.12e\n", ParticleNumber[n], Time,
    cell, ParticlePosition[0][n], ParticlePosition[1][n], ParticlePosition[2][n]);

//      printf("New star particle index = %d (%d %d)\n",
//	     ParticleNumber[n], NumberCount1, NumberCount2);
    }
  fclose(fp);
  return;
}
