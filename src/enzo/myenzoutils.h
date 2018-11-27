#ifndef MY_ENZO_UTILS_H
#define MY_ENZO_UTILS_H

#include "stddef.h"
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

template<class T>
inline const bool areOfSameType(T a, T b);

template<class T, class U>
inline const bool areOfSameType(T a, U b);

/**
 * Array allocation and initialization functions.
 */

template<typename T, typename U>
T* arr_cpy(T* dest, U* src, size_t n);

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

template<typename T, typename U>
T* arr_newset(size_t n, U x = 0);

template<typename T, typename U>
T* arr_newset(T** a, size_t n, U x);

template<typename T, typename U>
T* arr_set(T* a, size_t n, U x);

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
 * Math, BLAS-like and numpy-like functions.
 */
template<typename T>
inline T square(T x);

// Returns x[0] + ... + x[n-1]
template<typename T>
T arr_sum(T* x, const size_t n);

// x[i] := sum(x[0], ..., x[i])
template<typename T, typename U>
T* arr_cumsum(T* dest, size_t n, U initialSum);

// y[i] := sum(x[0], ..., x[i])
template<typename T, typename U, typename V>
T* arr_cumsum(T* dest, U* x, size_t n, V initialSum);

// y := x+y
// no overlapping
template<typename T, typename U>
T* arr_xpy(T* dest, U* x, const size_t n);

// x := a*x
template<typename T, typename U>
T* arr_ax(T* x, const size_t n, const U a);

// y := a*x + y
//no overlapping
template<typename T, typename U, typename V>
T* arr_axpy(T* dest, U* x, const size_t n, const V a);

// y := a*x + b*y
//no overlapping
template<typename T, typename U, typename V>
T* arr_axpby(T* dest, U* x, const size_t n, const V a);

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

#include "myenzoutils.T"

//#end if /*__macros_and_parameters_h_ */
#endif /* MY_ENZO_UTILS_H */
