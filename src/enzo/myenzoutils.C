#include "stddef.h"

#include "math.h"
#include "myenzoutils.h"

inline long double distancel(long double x1, long double x2, long double y1, long double y2)
{
	return sqrtl(square(y1 - x1) + square(y2 - x2));
}

inline long double distancel(long double x1, long double x2, long double x3, long double y1, long double y2,
								long double y3)
{
	return sqrtl(square(y1 - x1) + square(y2 - x2) + square(y3 - x3));
}

int mpiErrorString(char* s, int mpiError)
{
	int msglen;
	if(0 <= mpiError && mpiError <= MPI_ERR_LASTCODE)
		return MPI_Error_string(mpiError, s, &msglen);
	s[0] = '\0';
	return 0;
}

int mpiWait(MPI_Request* request)
{
	MPI_Status status;
	status.MPI_ERROR = NO_SUCH_MPI_ERROR;
	status.MPI_SOURCE = -1;
	status.MPI_TAG = -1;

	MPI_Wait(request, &status);

	return (status.MPI_ERROR == NO_SUCH_MPI_ERROR) ? 0 : status.MPI_ERROR;
}
