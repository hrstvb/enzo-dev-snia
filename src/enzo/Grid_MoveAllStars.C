/***********************************************************************
/
/  GRID CLASS (MOVE ALL STARS FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
/    NOTE: Adapted from grid::MoveAllParticles()
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void DeleteStarList(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);
Star *PopStar(Star * &Node);
Star* StarBufferToList(StarBuffer *buffer, int n);

int grid::MoveAllStars(int NumberOfGrids, grid* FromGrid[], int TopGridDimension)
{

  Star *NewStar, *cstar;

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }

  /* Determine total number of stars. */

  int grid_num;
  int TotalNumberOfStars = NumberOfStars;

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++)
    TotalNumberOfStars += FromGrid[grid_num]->NumberOfStars;
  if (TotalNumberOfStars == 0)
    return SUCCESS;

  /* Determine level of this grid */

  int ThisLevel = nint(-logf(TopGridDimension * CellWidth[0][0]) / M_LN2);

  /* Debugging info. */

//  if (debug) printf("MoveAllStars: %"ISYM" (before: ThisGrid = %"ISYM").\n",
//		    TotalNumberOfStars, NumberOfStars);

  /* Copy FromGrids' stars to new grid. */

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {

   /* Move to grid but not to the host processor.  We will take care
      of that in CommunicationCollectParticles. */

    if (FromGrid[grid_num]->Stars != NULL) {
      cstar = FromGrid[grid_num]->Stars;
      while (cstar != NULL) {
	NewStar = PopStar(cstar); // also advances cstar pointer
	NewStar->CurrentGrid = this;
	NewStar->level = ThisLevel;
	NewStar->GridID = this->ID;
	InsertStarAfter(this->Stars, NewStar);
	//cstar = cstar->NextStar;
      } // ENDWHILE stars
      FromGrid[grid_num]->Stars = NULL;
    }

//    /* Otherwise, communicate. */
//
//    else {
//      if (MyProcessorNumber == ProcessorNumber ||
//          MyProcessorNumber == FromGrid[grid_num]->ProcessorNumber)
//	if (FromGrid[grid_num]->CommunicationSendStars(this, ProcessorNumber) == FAIL) {
//	  ENZO_FAIL("Error in grid->CommunicationSendStars.\n");
//        }
//
//    } // ENDELSE same processor

  } // end: loop over grids.

  /* Set new number of stars in this grid. */

  NumberOfStars = TotalNumberOfStars;

  /* Set number of stars to zero.  No need to delete particles
     here, which is done in CommunicationSendStars.  For local copy,
     we copy the pointer to the new grid, and we shouldn't delete the
     linked list.  */

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    FromGrid[grid_num]->NumberOfStars = 0;
//    if (MyProcessorNumber == FromGrid[grid_num]->ProcessorNumber)
//      DeleteStarList(FromGrid[grid_num]->Stars);
  }

  return SUCCESS;
}












/* Old version of MoveAllStars : 
   might just work for less computationally intensive runs */

int grid::MoveAllStarsOld(int NumberOfGrids, grid* FromGrid[], int TopGridDimension)
{

  Star *NewStar, *cstar;

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }

  /* Determine total number of stars. */

  int grid_num;
  int TotalNumberOfStars = NumberOfStars;

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++)
    TotalNumberOfStars += FromGrid[grid_num]->NumberOfStars;
  if (TotalNumberOfStars == 0)
    return SUCCESS;

  /* Determine level of this grid */

  int ThisLevel = nint(-logf(TopGridDimension * CellWidth[0][0]) / M_LN2);

  /* Debugging info. */

//  if (debug) printf("MoveAllStars: %"ISYM" (before: ThisGrid = %"ISYM").\n",
//		    TotalNumberOfStars, NumberOfStars);

  /* Copy FromGrids' stars to new grid. */

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {

   /* Move to grid but not to the host processor.  We will take care
      of that in CommunicationCollectParticles. */

    if (FromGrid[grid_num]->Stars != NULL) {
      cstar = FromGrid[grid_num]->Stars;
      while (cstar != NULL) {
	NewStar = PopStar(cstar); // also advances cstar pointer
	NewStar->CurrentGrid = this;
	NewStar->level = ThisLevel;
	NewStar->GridID = this->ID;
	InsertStarAfter(this->Stars, NewStar);
	//cstar = cstar->NextStar;
      } // ENDWHILE stars
      FromGrid[grid_num]->Stars = NULL;
    }

    /* Otherwise, communicate. */

    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[grid_num]->ProcessorNumber)
	if (FromGrid[grid_num]->CommunicationSendStars(this, ProcessorNumber) == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendStars.\n");
        }

    } // ENDELSE same processor

  } // end: loop over grids.

  /* Set new number of stars in this grid. */

  NumberOfStars = TotalNumberOfStars;

  /* Set number of stars to zero.  No need to delete particles
     here, which is done in CommunicationSendStars.  For local copy,
     we copy the pointer to the new grid, and we shouldn't delete the
     linked list.  */

  for (grid_num = 0; grid_num < NumberOfGrids; grid_num++) {
    FromGrid[grid_num]->NumberOfStars = 0;
    if (MyProcessorNumber == FromGrid[grid_num]->ProcessorNumber)

      DeleteStarList(FromGrid[grid_num]->Stars);
  }

  return SUCCESS;
}
