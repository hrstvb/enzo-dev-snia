#include "stddef.h"

#include "IDE_defs.h"
#include "math.h"
#include "myenzoutils.h"

#ifdef USE_MPI

/*
 * Returns the MPI error string if defined for this mpiError,
 * otherwise returns the empty string.
 */
int mpiErrorString(char* s, int mpiError)
{
	int msglen;
	if(0 <= mpiError && mpiError <= MPI_ERR_LASTCODE)
		return MPI_Error_string(mpiError, s, &msglen);
	s[0] = '\0';
	return 0;
}

/*
 * An MPI_Wait wrapper returning the error code
 * from the received status.
 */
int mpiWait(MPI_Request* request)
{
	MPI_Status status;
	status.MPI_ERROR = NO_SUCH_MPI_ERROR;
	status.MPI_SOURCE = -1;
	status.MPI_TAG = -1;

	MPI_Wait(request, &status);

	return (status.MPI_ERROR == NO_SUCH_MPI_ERROR) ? 0 : status.MPI_ERROR;
}

#endif /* USE_MPI */
