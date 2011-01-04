/// @file      AMRsolve_domain.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Domain class

#include <stdio.h>
#include <assert.h>

#include <string>

#include "AMRsolve_scalar.h"
#include "AMRsolve_domain.h"

//----------------------------------------------------------------------

AMRsolve_Domain::AMRsolve_Domain() throw() : d_(0)
{ 
  for (int i=0; i<3; i++) xl_[i] = xu_[i] = 0.0;
}

//----------------------------------------------------------------------


AMRsolve_Domain::AMRsolve_Domain(std::string parms) throw() : d_(0)
{
  // Define a domain given text parameters, typically from an input file
  input (parms);
  assert (1 <= d_ && d_ <= 3);
}

//----------------------------------------------------------------------


AMRsolve_Domain::AMRsolve_Domain(int d, Scalar xl[3], Scalar xu[3]) 
  throw() : d_(d)
{
  assert (1 <= d_ && d_ <= 3);
  for (int i=0; i<d; i++) {
    xl_[i] = xl[i];
    xu_[i] = xu[i];
  }
}

//----------------------------------------------------------------------

AMRsolve_Domain::~AMRsolve_Domain() throw()
{
  //
}
  
//======================================================================

void AMRsolve_Domain::print() throw()
{
  printf("AMRsolve_Domain\n"
	 "   dimension      %d\n"
	 "   lower position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
	 "   upper position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n",
	 d_,
	 xl_[0],xl_[1],xl_[2],
	 xu_[0],xu_[1],xu_[2]);
}

//----------------------------------------------------------------------

void AMRsolve_Domain::write(FILE *fp) throw()
{
  assert (d_ == 3);
  if (fp == 0) fp = stdout;

  fprintf(fp,"domain "
	  "%d "
	  SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF
	  SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF,
	  d_,
	  xl_[0],xl_[1],xl_[2],
	  xu_[0],xu_[1],xu_[2]);
}

//----------------------------------------------------------------------

void AMRsolve_Domain::input(std::string parms) throw()
{
  sscanf(parms.c_str(),
	 "%d"
	 SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF
	 SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF,
	 &d_,&xl_[0],&xl_[1],&xl_[2],&xu_[0],&xu_[1],&xu_[2]);
}

//----------------------------------------------------------------------
