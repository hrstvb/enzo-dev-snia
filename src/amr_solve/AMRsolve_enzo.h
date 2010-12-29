// This file is for including Enzo grid class declarations in AMRsolve* files
// jobordner@ucsd.edu
// Updated: Daniel Reynolds (reynolds@smu.edu)

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "ExternalBoundary.h"
#include "ProtoSubgrid.h"
#include "GridList.h"
#include "Grid.h"
#include "LevelHierarchy.h"
#include "Hierarchy.h"

#ifdef int
#  undef int
#endif

#ifdef float
#  undef float
#endif
