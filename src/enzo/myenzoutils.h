#ifndef MY_ENZO_UTILS_H
#define MY_ENZO_UTILS_H

#include "stddef.h"
#include "math.h"
#include "string.h"
#include "IDE_defs.h"
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

using namespace std;

#ifdef __cplusplus
#define STD_CPP __cplusplus
#endif /* __cplusplus */

#define STD_CPP98 199711L
#define STD_CPP11 201103L
#define STD_CPP14 201402L
#define STD_CPP17 201703L

#ifdef __STDC__
#ifdef __STDC_VERSION__
#define STD_C __STDC_VERSION__
#else
#define STD_C
#endif /* __STDC_VERSION__ */
#endif /* __STDC__ */

#define STD_C95 199409L
#define STD_C99 199901L
#define STD_C11 201112L
#define STD_C17 201710L

#ifndef M_PI
# define M_PI		3.14159265358979323846
#endif
#ifndef M_PIl
# define M_PIl		3.141592653589793238462643383279502884L
#endif
#ifndef M_GOLDEN_RATIO
#define M_GOLDEN_RATIO      1.61803398874989
#endif
/*
 * Type utilities.
 */
template<class T>
inline const bool areOfSameType(T a, T b);

template<class T, class U>
inline const bool areOfSameType(T a, U b);

/**
 * Array allocation and initialization function templates.
 */

template<typename T>
inline T* arr_delbrnew(T** a, size_t n);

template<typename T, typename U>
T* arr_delbrnewcpy(T** a, U* xarr, size_t n);

template<typename T, typename U>
T* arr_delbrnewset(T** a, size_t n, U x);

template<typename T>
inline T* arr_delnew(T** a, size_t n);

template<typename T, typename U>
T* arr_delnewcpy(T** a, U* xarr, size_t n);

template<typename T, typename U>
T* arr_delnewset(T** a, size_t n, U x);

template<typename T, typename U>
T* arr_newcpy(U* xarr, size_t n);

template<typename T, typename U>
T* arr_newcpy(T** a, U* xarr, size_t n);

//template<typename T, typename U>
//T* arr_newset(size_t n, U x);
template<typename T>
T* arr_newset(size_t n, T x);

template<typename T, typename U>
T* arr_newset(T** a, size_t n, U x);

// Cannot overlap, if T != U.
template<typename T, typename U>
T* arr_cpy(T* dest, U* src, size_t n);

template<typename T, typename U>
T* arr_set(T* a, size_t n, U x);

template<typename T>
T** arr_set(T** a, size_t n, void* x);

template<typename T>
T** arr_set(T** a, size_t n, signed char x);
template<typename T>
T** arr_set(T** a, size_t n, short x);
template<typename T>
T** arr_set(T** a, size_t n, int x);
template<typename T>
T** arr_set(T** a, size_t n, long x);
template<typename T>
T** arr_set(T** a, size_t n, long long x);
template<typename T>
T** arr_set(T** a, size_t n, unsigned char x);
template<typename T>
T** arr_set(T** a, size_t n, unsigned short x);
template<typename T>
T** arr_set(T** a, size_t n, unsigned int x);
template<typename T>
T** arr_set(T** a, size_t n, unsigned long x);
template<typename T>
T** arr_set(T** a, size_t n, unsigned long long x);

//typedef int new_delop_t;
//const new_delop_t new_no_del = 0;
//const new_delop_t new_del = 1;
//const new_delop_t new_del_brackets = 2;
//
/////*
// * Redefined new operator with the same functionality.
// */
//template<typename T>
//void *operator new(size_t size, T** a, new_delop_t delete_op);
//
//template<typename T>
//void *operator new(size_t size, T** a, new_delop_t delete_op, T x);

/*
 * array search
 */
template<class T, class U>
size_t findmaxlte(T* a, size_t n, T x);

/*
 * Math, BLAS-like and numpy-like function templates.
 */

template<typename T, typename U>
inline long double lenl(T* x[], U* y[], size_t n);

template<typename T>
inline T square(T x);

template<typename T>
inline T cube(T x);

inline long double lenl(long double x1, long double x2, long double y1, long double y2);

inline long double lenl(long double x1, long double x2, long double x3, long double y1, long double y2,
								long double y3);

inline long double lenl(long double x1, long double x2, long double x3);

template<typename T, typename U>
long double lenl(T* x, U* y, size_t n);

template<typename T, typename U>
long double lenl(T* x, size_t n);

inline long double lensqaredl(long double x1, long double x2, long double y1, long double y2);

inline long double lensqaredl(long double x1, long double x2, long double x3, long double y1, long double y2,
								long double y3);
inline long double lensqaredl(long double x1, long double x2, long double x3);

template<typename T, typename U>
long double lensqaredl(T* x, U* y, size_t n);

template<typename T, typename U>
long double lensqaredl(T* x, size_t n);

template<typename T, typename U>
int normalizel(T* x, size_t n);

template<typename T>
inline int sign(T x)
{
	return (x > 0) - (x < 0);
}

// Returns x[0] + ... +, x[n-1]
template<typename T>
T arr_sum(T* x, const size_t n);

// x[i] := cumsum(x[0], ..., x[i])
template<typename T, typename U>
T* arr_cumsum(T* x, size_t n, U initialSum);

// y[i] := cumsum(x[0], ..., x[i])
// no overlapping
template<typename T, typename U, typename V>
T* arr_cumsum(T* dest, U* x, size_t n, V initialSum);

// y := x+y
// no overlapping if T != U
template<typename T, typename U>
T* arr_xpy(T* dest, const U* x, const size_t n);

// x := a*x
template<typename T, typename U>
T* arr_ax(T* x, const size_t n, const U a);

// dest := a*x
template<typename T, typename U, typename V>
T* arr_ax(T* dest, const U* x, const size_t n, const V a);

// y := a*x + y
//no overlapping if T != U
template<typename T, typename U, typename V>
T* arr_axpy(T* dest, const U* x, const size_t n, const V a);

// y := a*x + b*y
//no overlapping
template<typename T, typename U, typename V, typename W>
T* arr_axpby(T* dest, U* x, const size_t n, const V a, const W b);

template<typename T>
size_t sprintfvec(char* const s, const char* const preffix, const char* const singleElementFormat,
	const char* const separator, const char* const suffix, const T* const vec, const size_t n,
	const bool printWithIndices, const bool indexPrintsBeforeElement);

/*
 * Segmented array routines
 *
 * A segmented array consists of alternating segments of
 * data and undefined values.
 */
template<typename T>
void segarrPackByDefaultValue(T* packedarr, size_t* seglens, size_t* packedcount, size_t* segcount, size_t n, T* a,
								T defaultValue = (T) 0);

size_t segarrcounts(size_t* seglens, size_t* segcount);

template<typename T>
size_t segarrunpack(T* dest, T* packedarr, size_t* seglens, T defaultValue = (T) 0);

/*
 * MPI related
 */
#ifdef USE_MPI
#define NO_SUCH_MPI_ERROR (-1)

/*
 * mpiErrorString returns an empty string when the error code
 * would cause an error in MPI_Error_string.
 */
int mpiErrorString(char* s, int mpiError);

/*
 * Returns a non-zero error code of the request status
 * or 0 if no error.
 */
int mpiWait(MPI_Request* request);

/*
 * Returns an MPI_Datatype that corresponds to T, e.g.
 * @getMPI_Datatype<double>() -> MPI_DOUBLE.
 */
template<typename T>
inline MPI_Datatype getMPI_Datatype();

/*
 * Returns an MPI_Datatype that corresponds to T, e.g.
 * @double x;
 * @getMPI_Datatype<double>(x) -> MPI_DOUBLE.
 */
template<typename T>
inline MPI_Datatype getMPI_Datatype(T a);

/*
 * Returns an MPI_Datatype that corresponds to T, e.g.
 * @double* x;
 * @getMPI_Datatype<double>(x) -> MPI_DOUBLE.
 */
template<typename T>
inline MPI_Datatype getMPI_Datatype(T* a);
#endif /* USE_MPI */

int snlprintf(char* const s, const size_t size, size_t* const length, const char* const format, ...);

#include "myenzoutils.T"

#endif /* MY_ENZO_UTILS_H */
