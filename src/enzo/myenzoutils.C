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

/*
 * Similar to snprintf, but intended for successive appendd to
 * the same buffer.  The first call should be with length=0;
 *   const size_t N = 4;
 *   char s[N];
 *   int n, len = 0;
 *   n = snlprintf(s, N, &len, "abc");
 *   // n <- 3, len <- 3
 *
 * Subsequent call attempt to write at position s+length
 * and increase length with the number of number of characters
 * that were attempted to be written not including the
 * terminating null character.  The increment is returned
 * as the function result.  A negative value indicates an error.
 *
 *   n = snlprintf(s, N, &len, "def");
 *   // n <- 3, len <-6
 *   n = snlprintf(s, N, &len, "ghih");
 *   // n <- 3, len <-10
 *
 */
int snlprintf(char* const s, const size_t size, size_t* const length, const char* const format, ...)
{
	va_list varargs;
	int n;
	const size_t len = *length;

	va_start(varargs, format);
	if(size > len + 1)
		n = vsnprintf(s + len, size - len, format, varargs);
	else
		n = vsnprintf(NULL, 0, format, varargs);
	va_end(varargs);

	if(n > 0)
		*length += n;

	return n;
}


#endif /* USE_MPI */
