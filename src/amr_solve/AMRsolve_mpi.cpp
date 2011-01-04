/// @file      AMRsolve_mpi.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Mpi class

#include <stdio.h>
#include <mpi.h>
#include "AMRsolve_mpi.h"

AMRsolve_Mpi* pmpi = NULL;
