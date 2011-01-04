#ifdef AMR_SOLVE
#ifndef AMR_SOLVE_H
#define AMR_SOLVE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "HYPRE_sstruct_ls.h"

#include "AMRsolve_scalar.h"
#include "AMRsolve_point.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_error.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_hypre_grav.h"
#include "AMRsolve_hypre_fld.h"

#endif /* AMR_SOLVE_H (prevent multiple includes) */
#endif /* AMR_SOLVE (configuration) */
