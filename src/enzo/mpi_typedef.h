#ifndef __MPI_TYPEDEF_H
#define __MPI_TYPEDEF_H

typedef struct {
  char CallType;
  int GridNum[2];
  int GridOffset[3];
  int GridDims[3];
  FLOAT FArg[3];
  int IArg[3];
} mpi_header;

typedef struct {
  mpi_header header;
  void *buffer;
} enzo_message;

#endif /* ifndef __MPI_TYPEDEF_H */
