/*
 * IDE_defs.h
 *
 *  Created on: Nov 6, 2018
 *      Author: B. Hristov
 *
 *  Preprocessor definitions, which are normally defined on
 *  the compiler command-line.  It is intended for use in
 *  integrated development environments (IDEs) to enable
 *  proper syntax highlighting and functioning.  In order
 *  for these definitions to take effect one needs to define
 *  (in an IDE-specific way) the symbol IDE_PRESENT with a
 *  value of 1.
 *
 *  This file is not automatically generated.
 *
 *  See also: Makefile, $(DEFINES) for usage;
 *            auto_show_flags.C, "DEFINES = -D... -D..." for
 *            actual symbols and values.
 */
#ifndef IDE_DEFS_H
#define IDE_DEFS_H

#ifdef IDE_PRESENT
#if IDE_PRESENT == 1

#define LINUX
#define H5_USE_16_API
#define MAKING -w
#define __max_subgrids 100000
#define __max_baryons 30
#define __max_cpu_per_node 8
#define __memory_pool_size 100000
#define INITS64
#define LARGE_INTS
#define CONFIG_PINT_8
#define IO_32
#define USE_MPI
#define CONFIG_PFLOAT_8
#define CONFIG_BFLOAT_8
#define USE_HDF5_GROUPS
#define TRANSFER
#define NEW_GRID_IO
#define FAST_SIB
#define ENZO_PERFORMANCE
#define SAB

#endif /* IDE_PRESENT==1 */
#endif /* IDE_PRESENT */

#include "hdf5.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#endif /* IDE_DEFS_H */
