// $Id: newgrav-hypre-solve.h 10 2010-03-18 22:35:19Z jbordner $

/// @file      newgrav-hypre-solve.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Useful definitions for the hypre-solve project

#ifndef NEWGRAV_HYPRE_SOLVE_H
#define NEWGRAV_HYPRE_SOLVE_H

#define _TRACE_ if (trace) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TRACE %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }

#define _TRACE_HYPRE_ if (trace_hypre) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TRACE HYPRE%d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }
#define _BARRIER_ if (trace) { MPI_Barrier (MPI_COMM_WORLD); }

#define _TEMPORARY_ { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TEMPORARY %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
#define MAX(a,b) ( (a) > (b) ? (a) : (b))

#endif /* NEWGRAV_HYPRE_SOLVE_H */


