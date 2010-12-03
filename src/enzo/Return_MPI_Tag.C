#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

// Return a somewhat-unique MPI tag for communication.

MPI_Arg Return_MPI_Tag(int tag, int num1[], int num2[])
{
  //return 37*tag + 17*num1;// + 7*num2;
  //return 1000000*tag + num1;// + 7*num2;
  long_int result = 0;

  // Flag n-th bit (up to 64) for the tag type
  result |= (1 << (tag % 64));

  // Increase tag by the bits num1[0]
  result |= num1[0];

  // Increase tag by the bits num1[1], shifted by 16 bits
  result |= (num1[1] << 16);

  // If specified, increase by num2[], shifted by 8, 16, and 24 bits.
  // There will be some overlap, but the collision probability between
  // tags is low.
  if (num2 != 0) {
    result |= (num2[0] << 8);
    result |= (num2[1] << 16);
    result |= (num2[2] << 24);
  }

  // Reduce to a 32-bit signed integer
  result = result % (1 << 31);

  if ((MPI_Arg) result < 0) {
    printf("MPI_Tag < 0 !!! :: result = %d, tag = %d, num1 = %d %d",
	   result, tag, num1[0], num1[1]);
    if (num2 != 0)
      printf(", num2 = %d %d %d", num2[0], num2[1], num2[2]);
    printf("\n");
    ENZO_FAIL("");
  }

  return (MPI_Arg) result;
  
}
