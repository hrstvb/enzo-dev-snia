/// @file      newgrav-mpi.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Mpi class

#ifndef AMRSOLVE_MPI_H
#define AMRSOLVE_MPI_H

#ifndef TRACE
#   define TRACE {printf("TRACE %s:%d\n",__FILE__,__LINE__); fflush(stdout);}
#endif

class AMRsolve_Mpi {

  /// @class    AMRsolve_Mpi
  /// @brief    AMRsolve_Mpi class of MPI convenience functions

public:

  AMRsolve_Mpi() : comm_(0), np_(1), ip_(0) 
  { }

  AMRsolve_Mpi(MPI_Comm comm)
  { comm_ = comm; initialize_(); }

  AMRsolve_Mpi(int* argc, char*** argv) : comm_(MPI_COMM_WORLD)
  {
    TRACE;
    printf("argc=%d  argv[0]=%s\n",*argc, *argv[0]);
    MPI_Init(argc,argv);
    TRACE;
    initialize_();
    TRACE;
  }

  ~AMRsolve_Mpi()
  {
    //    MPI_Finalize (); // MPI_Finalize() seems to be called at program exit
    //                        complains if included
  };

  bool is_root() throw() {return ip_ == 0;};
  int ip() throw() {return ip_;};
  int np() throw() {return np_;};
  void barrier() throw() { MPI_Barrier(comm_); };
  void initialize_()
  { MPI_Comm_size(comm_, &np_); MPI_Comm_rank (comm_, &ip_); }
    
private:

  MPI_Comm comm_;  // communicator
  int np_;         // Number of processors
  int ip_;         // Rank of this processor

};

extern AMRsolve_Mpi * pmpi;

#endif /* AMRSOLVE_MPI_H */
