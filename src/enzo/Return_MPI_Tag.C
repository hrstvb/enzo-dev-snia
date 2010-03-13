#include <stdio.h>
#include "macros_and_parameters.h"

// Return a somewhat-unique MPI tag for communication.  The factors
// are prime.

MPI_Arg Return_MPI_Tag(int tag, int num1, int num2)
{
  return 37*tag + 17*num1;// + 7*num2;
}
