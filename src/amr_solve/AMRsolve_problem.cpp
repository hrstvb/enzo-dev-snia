/// @file      problem.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Problem class

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <map>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "AMRsolve_defs.h"


//----------------------------------------------------------------------

#define debug false
#define trace false

//----------------------------------------------------------------------

#include "AMRsolve_scalar.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_error.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_point.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_problem.h"

//======================================================================

AMRsolve_Problem::AMRsolve_Problem() throw()
{
  //
}

//----------------------------------------------------------------------

AMRsolve_Problem::~AMRsolve_Problem() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

AMRsolve_Problem::AMRsolve_Problem(const AMRsolve_Problem& p) throw()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
}

//----------------------------------------------------------------------

AMRsolve_Problem& AMRsolve_Problem::operator=(const AMRsolve_Problem& p) throw()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
  return *this;
}

//----------------------------------------------------------------------

void AMRsolve_Problem::read(std::string filename) throw()
{
  bool boundary_defined = false;

  parameters_.read(filename);
  ItParameters itp (parameters_);

  while (itp++) {

    std::string key   = itp.key();
    std::string value = itp.value();

    if (key == "dimension") {
      int d = atoi(value.c_str());
      hierarchy_.set_dim(d);
      AMRsolve_Point::set_dim(d);

    // Domain
    } else if (key == "domain") {
      domain_.input(value);

    // Grid ...
    } else if (key == "grid") {
      hierarchy_.insert_grid(new AMRsolve_Grid(value));

    } else if (key == "point") {
      points_.push_back(new AMRsolve_Point(value));      

    } else if (key == "boundary") {
      boundary_defined = true;
    }
  }

  if (! boundary_defined ) {
    fprintf (stderr,"Input parameter 'boundary' must be defined\n");
    MPI_Finalize();
    exit(1);
  }
}

//----------------------------------------------------------------------

void AMRsolve_Problem::print() throw()
{
  int i;
  domain_.print();
  hierarchy_.print();
  for (i=0; i<num_points(); i++)  point(i).print();
}

//----------------------------------------------------------------------

void AMRsolve_Problem::write(FILE *fp) throw()
{
  if (fp == 0) fp = stdout;
  int i;
  domain_.write(fp);
  hierarchy_.write(fp);
  for (i=0; i<num_points(); i++)  point(i).write(fp);
}

//----------------------------------------------------------------------

void AMRsolve_Problem::deallocate_() throw()
{
  for (unsigned i=0; i<points_.size(); i++) {
    delete points_[i];
    points_[i] = 0;
  }
  points_.resize(0);
}

//----------------------------------------------------------------------
