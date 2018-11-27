#include "stddef.h"

#include "myenzoutils.h"

/*
 * Array (re)allocation fuctions @arr...() with
 * initialization.  All functions return a reference to the
 * new array or NULL if @new coudln't allocate the requested
 * memeory.
 *
 * The @arrdel...() functions perform @delete @*a in the
 * begining, unless @*a==NULL.  Similarly the functions
 * named @arrdelbrnew...() use the delete[] operator.
 *
 * The requested number of elements is @n. The new array is
 * allocated with @new.  If @a contains a reference to a
 * pointer, the latter is updated to point to the new array.
 * If @a == NULL, it is ignored.
 *
 * The @arr...set() functions set each element of the new
 * array to @x upon successsful allocation.
 *
 * Errors and exceptions may come from the delete operation,
 * for example if an object had been deleted already or if
 * a pointer is invalid.
 */

#define NO_SUCH_MPI_ERROR (-1)
int mpiErrorString(char* s, int mpiError)
{
	int msglen;
	if (0 <= mpiError && mpiError <= MPI_ERR_LASTCODE)
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
