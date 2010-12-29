/// @file      level.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Level class

#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <string>
#include <vector>

#include <mpi.h>

//----------------------------------------------------------------------

const int debug = 0;

#include "AMRsolve_defs.h"

#include "AMRsolve_scalar.h"
#include "AMRsolve_error.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"

//----------------------------------------------------------------------

AMRsolve_Domain AMRsolve_Level::domain_;

//----------------------------------------------------------------------

AMRsolve_Level::AMRsolve_Level(int n) throw() : n_(n)
{
  // Initialize extents
  for (int i=0; i<3; i++) {
    il_[i] = INT_MAX;
    iu_[i] = INT_MIN;
  }
  // Initialize grid list with 0 sentinal
  grids0_.push_back (0);
}
	  
//----------------------------------------------------------------------

AMRsolve_Level::~AMRsolve_Level() throw()
{
  // Grid objects deleted in hierarchy
  // Assumes all grids in all levels are in a unique hierarchy
  grids0_.resize (0);
  grids0_.push_back (0);
}

//======================================================================

void AMRsolve_Level::insert_grid(AMRsolve_Grid& grid) throw()
{
  // Append Grid to list of all Grids in the AMRsolve_Level
  grids0_[grids0_.size()-1] = & grid;
  grids0_.push_back (0);

  // Update AMRsolve_Level extents given Grid indices
  for (int i=0; i<3; i++) {
    il_[i] = MIN(il_[i],grid.index_lower(i));
    iu_[i] = MAX(iu_[i],grid.index_upper(i));
  }
}

//======================================================================

void AMRsolve_Level::print() throw()
{
  printf("AMRsolve_Level %d\n",n_);
  printf("   grids       = %d\n",int(grids0_.size()-1));
  printf("   index_lower = (%d,%d,%d)\n",il_[0],il_[1],il_[2]);
  printf("   index_upper = (%d,%d,%d)\n",iu_[0],iu_[1],iu_[2]);
  for (int i=0; i<num_grids(); i++) {
    return_grid(i).print();
  }
}

//----------------------------------------------------------------------

/// Write the AMRsolve_Level grids to the given open file in geomview format
void AMRsolve_Level::geomview_grid(FILE* fpr, bool full) throw()
{
  int n = num_grids();

  if (debug) printf("DEBUG %s:%d  geomview_grid n = %d\n",__FILE__,__LINE__,n);

  // Write AMRsolve_Level header
  if (full) {
    fprintf(fpr,"VECT\n");
    fprintf(fpr,"%d %d 1\n",2+4*n, 2+16*n);

    // Write vertices per polygon
    fprintf(fpr,"1 1 ");
    for (int i=0; i<n; i++) {
      fprintf(fpr,"8 3 3 2 ");
    }
    fprintf(fpr,"\n");

    // Write colors changes
    fprintf(fpr,"1 0 ");
    for (int i=0; i<n; i++) fprintf(fpr,"0 0 0 0 "); fprintf(fpr,"\n");

    // Print points at domain boundaries to provide geomview with bounding box

    Scalar dl[3],du[3];
    AMRsolve_Level::domain_.lower(dl[0],dl[1],dl[2]);
    AMRsolve_Level::domain_.upper(du[0],du[1],du[2]);
    fprintf(fpr,"%g %g %g\n",dl[0],dl[1],dl[2]);
    fprintf(fpr,"%g %g %g\n",du[0],du[1],du[2]);
  }

  // Write grids
  ItLevelGridsAll itg (*this);
  while (AMRsolve_Grid * g = itg++) {
    g->geomview_grid(fpr,0);
  }

  // Write trailer
  if (full) {
    fprintf(fpr,"1 1 1 0\n");
  }
}

//----------------------------------------------------------------------

/// Write the local AMRsolve_Level grids to the given open file in geomview format
void AMRsolve_Level::geomview_grid_local(FILE* fpr, bool full) throw()
{
  // Count local grids
  int n=0;
  ItLevelGridsLocal itg (*this);
  while (itg++) n++;

  if (debug) printf("DEBUG %s:%d  geomview_grid_local n = %d\n",__FILE__,__LINE__,n);

  // Write AMRsolve_Level header
  if (full) {
    fprintf(fpr,"VECT\n");
    fprintf(fpr,"%d %d 1\n",2+4*n, 2+16*n);

    // Write vertices per polygon
    fprintf(fpr,"1 1 ");
    for (int i=0; i<n; i++) {
      fprintf(fpr,"8 3 3 2 ");
    }
    fprintf(fpr,"\n");

    // Write colors changes
    fprintf(fpr,"1 0 ");
    for (int i=0; i<n; i++) fprintf(fpr,"0 0 0 0 "); fprintf(fpr,"\n");

    // Print points at domain boundaries to provide geomview with bounding box
    Scalar dl[3],du[3];
    AMRsolve_Level::domain_.lower(dl[0],dl[1],dl[2]);
    AMRsolve_Level::domain_.upper(du[0],du[1],du[2]);
    fprintf(fpr,"%g %g %g\n",dl[0],dl[1],dl[2]);
    fprintf(fpr,"%g %g %g\n",du[0],du[1],du[2]);
  }

  // Write grids
  while (AMRsolve_Grid * g = itg++) {
    g->geomview_grid (fpr,0);
  }

  // Write trailer
  if (full) {
    fprintf(fpr,"1 1 1 0\n");
  }
}

//----------------------------------------------------------------------

void AMRsolve_Level::geomview_face(FILE* fpr, bool full) throw()
{
  if (full) {
    fprintf(fpr,"CQUAD\n");
    // Print points at domain boundaries to provide geomview with bounding box

    Scalar dl[3],du[3];
    AMRsolve_Level::domain_.lower(dl[0],dl[1],dl[2]);
    AMRsolve_Level::domain_.upper(du[0],du[1],du[2]);
    fprintf(fpr,
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1\n",
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2]);
    fprintf(fpr,
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1\n",
	    du[0],du[1],du[2],
	    du[0],du[1],du[2],
	    du[0],du[1],du[2],
	    du[0],du[1],du[2]);
  }

  ItLevelGridsAll itga (*this);
  while (AMRsolve_Grid * grid = itga++) {
    if (debug) printf("DEBUG %s:%d grid id = %d\n",__FILE__,__LINE__,grid->id());
    grid->geomview_face(fpr,false);
  }
}

//----------------------------------------------------------------------

void AMRsolve_Level::geomview_face_types(FILE* fpr, 
					 AMRsolve_Faces::Label* types, 
					 int num_types, bool full) throw()
{
  if (full) {
    fprintf(fpr,"CQUAD\n");

    // Print points at domain boundaries to provide geomview with bounding box
    Scalar dl[3],du[3];
    AMRsolve_Level::domain_.lower(dl[0],dl[1],dl[2]);
    AMRsolve_Level::domain_.upper(du[0],du[1],du[2]);
    fprintf(fpr,
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1\n",
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2],
	    dl[0],dl[1],dl[2]);
    fprintf(fpr,
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1 "
	    "%g %g %g 0 0 0 1\n",
	    du[0],du[1],du[2],
	    du[0],du[1],du[2],
	    du[0],du[1],du[2],
	    du[0],du[1],du[2]);
  }

  ItLevelGridsAll itga (*this);
  while (AMRsolve_Grid * grid = itga++) {
    if (debug) printf("DEBUG %s:%d grid id = %d\n",__FILE__,__LINE__,grid->id());
    grid->geomview_face_type(fpr,types,num_types,false);
  }
}

//----------------------------------------------------------------------

void AMRsolve_Level::write(FILE *fp) throw()
{
  if (fp == 0) fp = stdout;

  for (int i=0; i<num_grids(); i++) {
    fprintf(fp,"AMRsolve_Level %d\n",n_);
    return_grid(i).write("header",fp);
    return_grid(i).write("u",fp);
    return_grid(i).write("f",fp);
  }
}

//----------------------------------------------------------------------
