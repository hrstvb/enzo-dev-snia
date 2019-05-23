#ifndef MY_ENZO_UTIILS_TPP
#define MY_ENZO_UTIILS_TPP
#include "stddef.h"
#include "math.h"
#include <stdexcept>
#include "IDE_defs.h"
#include "myenzoutils.h"
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

template<class T>
inline const bool areOfSameType(T a, T b)
{
	return true;
}
template<class T, class U>
inline const bool areOfSameType(T a, U b)
{
	return false;
}

/*
 * Array (re)allocation fuctions @arr...() with
 * initialization.  All functions return a reference to the
 * new array or NULL if @new coudln't allocate the requested
 * memory.
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
 * Errors and exceptions may come from the delete or new
 * operations.
 */

template<typename T>
inline T* arr_delbrnew(T** a, size_t n)
{
	delete[] *a;
	return *a = new T[n];
}

template<typename T, typename U>
T* arr_delbrnewcpy(T** a, U* xarr, size_t n)
{
	return arr_cpy(arr_delbrnew(a, n), n, xarr);
}

template<typename T, typename U>
T* arr_delbrnewset(T** a, size_t n, U x)
{
	return arr_set(arr_delbrnew(a, n), n, (T) x);
}

template<typename T>
inline T* arr_delnew(T** a, size_t n)
{
	delete *a;
	return *a = new T[n];
}

template<typename T, typename U>
T* arr_delnewcpy(T** a, U* xarr, size_t n)
{
	return arr_cpy(arr_delnew(a, n), n, xarr);
}

template<typename T, typename U>
T* arr_delnewset(T** a, size_t n, U x)
{
	return arr_set(arr_delnew(a, n), n, (T) x);
}

template<typename T, typename U>
T* arr_newcpy(U* xarr, size_t n)
{
	return arr_cpy(new T[n], n, xarr);
}

template<typename T, typename U>
T* arr_newcpy(T** a, U* xarr, size_t n)
{
	return *a = arr_cpy(new T[n], n, xarr);
}

//template<typename T, typename U>
//T* arr_newset(size_t n, U x)
//{
//	return arr_set(new T[n], n, (T) x);
//}
//
template<typename T>
T* arr_newset(size_t n, T x)
{
	return arr_set(new T[n], n, x);
}

template<typename T, typename U>
T* arr_newset(T** a, size_t n, U x)
{
	return *a = arr_set(new T[n], n, (T) x);
}

template<typename T, typename U>
T* arr_cpy(T* dest, U* src, size_t n)
{
	if(areOfSameType(dest, src))
		if((void*) src == (void*) dest)
			return dest;
		else
			return (T*) memmove(dest, src, n * sizeof(T));

	U* src_n = src + n;
	while(src < src_n)
		*dest++ = (*src++);
	return dest - n;
}

template<typename T, typename U, typename V>
T* arr_cpya(T* dest, U* src, size_t n, V a)
{
	if(a == 0)
		return arr_set(dest, n, 0);

	if(areOfSameType(dest, src))
		if((void*) src == (void*) dest)
			return dest;
		else
		{
			memmove(dest, src, n * sizeof(T));
			if(a == 1)
				return dest;
			return arr_ax(dest, n, a);
		}

	U* src_n = src + n;
	while(src < src_n)
		*dest++ = a * (*src++);
	return dest - n;
}

template<typename T, typename U>
T* arr_set(T* a, size_t n, U x)
{
	static int firstZero = 1;
	static unsigned char nonBytewiseZero;

	if(sizeof(T) == 1)
		return (T*) memset(a, (int) (long long) x, n);

	if(x == (T) 0)
	{
		if(firstZero)
		{
			firstZero = 0;

			for(unsigned char* c = (unsigned char*) &x; c < (unsigned char*) (&x + 1); c++)
				if(nonBytewiseZero = *c)
					break;
		}

		if(!nonBytewiseZero)
			return (T*) memset(a, 0, n * sizeof(T));
	}

	T* a_i = a + n;
	while(a < a_i)
		*--a_i = (T) x;
	return a;
}

template<typename T>
T** arr_set(T** a, size_t n, void* x)
{
	if(x)
		throw invalid_argument("arr_set: Only null can be used as a void pointer value.");
	return (T**) memset(a, 0, n * sizeof(T*));
}

//#if __cplusplus < STD_CPP11 /* C++11 */
//// long long should work most places but can be improved.
//#define INTPTR_T long long
//#define UINTPTR_T unsigned long long
//
//#else
//
//#define INTPTR_T intptr_t
//#define UINTPTR_T uintptr_t
//
//#endif /* __cplusplus < 201103L /* C++11) */

template<typename T>
T** arr_set(T** a, size_t n, signed char x)
{
	return arr_set(a, n, (long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, short x)
{
	return arr_set(a, n, (long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, int x)
{
	return arr_set(a, n, (long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, long x)
{
	return arr_set(a, n, (long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, long long x)
{
	if(x)
		throw invalid_argument("arr_set: Only 0 can be used as an implicit integer pointer value.");
	return (T**) memset(a, 0, n * sizeof(T*));
}

template<typename T>
T** arr_set(T** a, size_t n, unsigned char x)
{
	return arr_set(a, n, (unsigned long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, unsigned short x)
{
	return arr_set(a, n, (unsigned long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, unsigned int x)
{
	return arr_set(a, n, (unsigned long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, unsigned long x)
{
	return arr_set(a, n, (unsigned long long) x);
}

template<typename T>
T** arr_set(T** a, size_t n, unsigned long long x)
{
	if(x)
		throw invalid_argument("arr_set: Only 0 can be used as an implicit integer pointer value.");
	return arr_set(a, n, 0LL);
}

template<typename T, typename U>
T* arr_xy(T* dest, const U* y, const size_t n)
{
	T* x_n = dest + n;
	while(dest < x_n)
		*dest++ *= *y++;
	return dest;
}

template<typename T, typename U>
T* bits2arr01(T* const dest, size_t n, U i)
{
	arr_set(dest, n, 0);
	for(T* a = dest; a < dest + n && i; a++, i <<= 1)
		*a = i & 1;
}

template<typename T>
T* bits2arrMasks(T* const dest, const size_t n, const T i)
{
	arr_set(dest, n, 0);
	unsigned long long mask = 1;
	for(T* a = dest; a < dest + n; a++)
	{
		*a = i & mask;
		mask >>= 1;
	}
}

//#undef INTPTR_T
//#undef UINTPTR_T

/*
 * New @new operators encapsulating the functionality of
 * the @arr... functions.  Problably these shouldn't be
 * used to avoid confusion of those who read the code
 * later.
 */
//template<typename T>
//void *operator new(size_t size, T** a, new_delop_t* delete_op)
//{
//	if (a)
//	{
//		switch (&delete_op)
//		{
//		case new_del:
//			delete *a;
//			break;
//		case new_del_brackets:
//			delete[] *a;
//			break;
//		default: //no delete
//		}
//	}
//
//	T * b = new (nothrow) T[size];
//	if (a) *a = b;
//	return b;
//}
//
//template<typename T>
//void *operator new(size_t size, T** a, new_delop_t* delete_op, T x)
//{
//	T * b = new (a, delete_op) T[size];
//	if (b) arr_set(b, size, x);
//	return b;
//}
/*
 * Return the index of the largest element in the array @a,
 * which is less than or equal to @x.
 * Return -1, if @x < @a[0].
 * @a should be strictly ascending.
 */
template<class T, class U>
size_t findmaxlte(T* a, size_t n, U x)
{
	n--;
	size_t l = 0, m = n / 2, r = n;
	T z;
	if(!areOfSameType(a, x))
		z = x;
	T& y = (areOfSameType(a, x)) ? x : z;
	while(l < m)
	{
		if(a[m] <= y)
			l = m;
		else
			r = m;

		m = (l + r) / 2;
	}
	if(l == 0)
		return -(y < a[0]);
	if(r == n)
		return r - (y < a[n]);
	return l;
}

template<class T, class U>
long long indexof(const T* const a, size_t n, U x)
{
	const T* b = a + n;
	while(a < b)
		if(*--b == x)
			return b - a;
	return -1;
}

/*
 * Math
 */
template<typename T, typename U>
inline long double lenl(T* x[], U* y[], size_t n)
{
	long double r_squared = 0;
	T* x_n = x + n;
	while(x < x_n)
		r_squared += square(*y++ - *x++);

	return sqrtl(r_squared);
}

template<typename T>
inline T square(T x)
{
	return x * x;
}

template<typename T>
inline T cube(T x)
{
	return x * x * x;
}

inline long double lenl(long double x1, long double x2, long double y1, long double y2)
{
	return sqrtl(lensqaredl(x1, x2, y1, y2));
}

inline long double lenl(long double x1, long double x2, long double x3, long double y1, long double y2,
	long double y3)
{
	return sqrtl(lensqaredl(x1, x2, x3, y1, y2, y3));
}

inline long double lenl(long double x1, long double x2, long double x3)
{
	return sqrtl(lensqaredl(x1, x2, x3));
}

inline long double lensqaredl(long double x1, long double x2, long double y1, long double y2)
{
	return square(y1 - x1) + square(y2 - x2);
}

inline long double lensqaredl(long double x1, long double x2, long double x3, long double y1, long double y2,
	long double y3)
{
	return square(y1 - x1) + square(y2 - x2) + square(y3 - x3);
}

inline long double lensqaredl(long double x1, long double x2, long double x3)
{
	return square(x1) + square(x2) + square(x3);
}

template<typename T, typename U>
long double lensqaredl(T* x, U* y, size_t n)
{
	long double d = 0;
	T* x_n = x + n;
	while(x < x_n)
		d += square(*x++ - *y++);
	return d;
}

template<typename T, typename U>
long double lenl(T* x, U* y, size_t n)
{
	return sqrt(lenl(x, y, n));
}

template<typename T, typename U>
int normalizel(T* x, size_t n)
{
	long double L = lensqaredl(x, n);
	if(L)
	{
		arr_ax(x, n, 1/L);
		return 1;
	}
	return 0;
}

template<typename T, typename U, typename V>
long double distMaxSquaredl(const T* x1, const U* x2, const V* y, const size_t n)
{
	long double d = 0;
	long double d1, d2;
	const V* y_n = y + n;
	while(y < y_n)
	{
		d1 = square(*x1++ - *y);
		d2 = square(*x2++ - *y++);
		d += (d2 > d1) ? d2 : d1;
	}
	return d;
}

template<typename T, typename U, typename V>
long double distMaxl(const T* x1, const U* x2, const V* y, const size_t n)
{
	return sqrt(distMaxSquaredl(x1, x2, y, n));
}

template<typename T, typename U>
bool intersectRectangles(T* ldest, T* rdest, U* lother, U* rother, size_t n)
{
	for(size_t dim = 0; dim < n; dim++)
		if(ldest[dim] > rother[dim])
			return true;
	for(size_t dim = 0; dim < n; dim++)
		if(rdest[dim] < lother[dim])
			return true;
	for(size_t dim = 0; dim < n; dim++)
		if(ldest[dim] < lother[dim])
			ldest[dim] = lother[dim];
	for(size_t dim = 0; dim < n; dim++)
		if(rdest[dim] > rother[dim])
			rdest[dim] = rother[dim];
	return false;
}

/*
 * BLAS-like and numpy-like
 */
template<typename T>
T arr_min(const T* x, const size_t n)
{
	T minx = *x;
	const T* const x_n = x + n;

	while(++x < x_n)
		if(minx > *x)
			minx = *x;

	return minx;
}

template<typename T>
T arr_max(const T* x, const size_t n)
{
	T maxx = *x;
	const T* const x_n = x + n;

	while(++x < x_n)
		if(maxx < *x)
			maxx = *x;

	return maxx;
}

template<typename T>
void arr_minmax(T* min_x, size_t* min_i, T* max_x, size_t* max_i, const T* const x, const size_t n)
{
	T miny, maxy;
	size_t j, minj, maxj;
	if(n > 0)
	{
		const T* const x_n = x + n;
		const T* y = x;
		minj = maxj = j = 0;
		miny = maxy = *y;

		while(++y < x_n)
		{
			j++;
			if(miny > *y)
			{
				minj = j;
				miny = *y;
			}
			if(maxy < *y)
			{
				maxj = j;
				maxy = *y;
			}
		}
	}
	else
	{
		minj = maxj = -1;
	}

	if(min_x)
		*min_x = miny;
	if(min_i)
		*min_i = minj;
	if(max_x)
		*max_x = maxy;
	if(max_i)
		*max_i = maxj;
}

// Returns x[0] + ... +, x[n-1]
template<typename T>
T arr_sum(T* x, const size_t n)
{
	T sum = 0;
	const T* const x_n = x + n;
	while(x < x_n)
		sum += *x++;
	return sum;
}

// x[i] := cumsum(x[0], ..., x[i])
template<typename T, typename U>
T* arr_cumsum(T* x, size_t n, U initialSum)
{
	const T* x_n = x + n;

	if(areOfSameType(x, initialSum) || initialSum)
	{
		T sum = initialSum;
		while(x < x_n)
			sum = (*x++ += sum);
	}
	else
		while(x < x_n)
			initialSum = (*x++ += initialSum);

	return x - n;
}

// y[i] := cumsum(x[0], ..., x[i])
// no overlapping
template<typename T, typename U, typename V>
T* arr_cumsum(T* dest, U* x, size_t n, V initialSum)
{
	const T* x_n = x + n;

	if(areOfSameType(dest, x))
	{
		if(x == dest)
			return arr_cumsum(dest, n, initialSum);

		if(x < dest && dest < x_n)
		{
			arr_cpy(dest, x, n);
			return arr_cumsum(dest, n, initialSum);
		}
	}

	if(areOfSameType(x, initialSum) || initialSum == 0)
	{
		U sum = initialSum;
		while(x < x_n)
			*dest++ = (sum += *x++);
		return dest - n;
	}

	while(x < x_n)
		*dest++ = (initialSum += *x++);
	return dest - n;
}

#define ARR_LINOPBODY(DEST_EQ_SRC, DEST_IN_SRC, ASC_LOOP) \
	U* x_n = x + n; \
	if (areOfSameType(dest, x)) \
	{ \
		if ((void*)x == (void*)dest) \
			return DEST_EQ_SRC; \
	\
		if ((void*)x < (void*)dest && (void*)dest < (void*)x_n) \
		{ \
			dest += n; \
			while (x < x_n) \
				DEST_IN_SRC; \
			return dest; \
		} \
	} \
	while (x < x_n) \
		ASC_LOOP; \
	return dest - n;

// y := x+y
// no overlapping if T != U
template<typename T, typename U>
T* arr_xpy(T* dest, U* x, const size_t n)
{
	ARR_LINOPBODY(arr_ax(dest, n, 2), *--dest += *--x_n, *dest++ += *x++)
}

// x := x + a
template<typename T, typename U>
T* arr_xpa(T* x, const size_t n, const U a)
{
	if(a == 0)
		return x;

	T* x_n = x + n;
	while(x < x_n)
		*--x_n += a;
	return x;
}

// x := a*x
template<typename T, typename U>
T* arr_ax(T* x, const size_t n, const U a)
{
	if(a == 1)
		return x;
	if(a == 0)
		return arr_set(x, n, (T) 0);
	T* x_n = x + n;
	while(x < x_n)
		*--x_n *= a;
	return x;
}

// x := a*x
template<typename T, typename U, typename V>
T* arr_axpb(T* x, const size_t n, const U a, const V b)
{
	if(a == 1)
		if(b==0)
			return x;
		else
			return arr_xpa(x, n, b);

	if(a == 0)
		if(b == 0)
			return arr_set(x, n, (T) 0);
		else
			return arr_set(x, n, (T) b);

	if(b==0)
		return arr_ax(x, n, a);

	T* x_n = x + n;
	while(x < x_n)
		*--x_n *= a;
	return x;
}

// dest := a*x
template<typename T, typename U, typename V>
T* arr_ax(T* dest, U* x, const size_t n, const V a)
{
	if(a == 0)
		return arr_set(dest, n, 0);
	if(a == 1)
		return arr_cpy(dest, x, n);

	ARR_LINOPBODY(arr_ax(x, n, a), *--dest = a * *--x, *dest++ = a * *x++)
}

// y := a*x + y
//no overlapping if T != U
template<typename T, typename U, typename V>
T* arr_axpy(T* dest, U* x, const size_t n, const V a)
{
	if(a == 0)
		return dest;
	if(a == 1)
		return arr_xpy(dest, x, n);

	ARR_LINOPBODY(arr_ax(dest, n, a + 1), *--dest += a * *--x_n, *dest++ += a * (*x++))
}

// y := a*x + b*y
template<typename T, typename U, typename V, typename W>
T* arr_axpby(T* dest, U* x, const size_t n, const V a, const W b)
{
	if(b == 0)
	{
		arr_cpy(dest, x, n);
		return arr_ax(dest, n, a);
	}
	if(a == 0)
		return arr_ax(dest, n, b);
	if(b == 1)
		return arr_axpy(dest, x, n, a);

	ARR_LINOPBODY(arr_ax(x, n, a + b),
	{
		*--dest *= b
		;
		*dest += a * (*--x_n)
		;
	}
					,
					{
						*dest *= b
						;
						*dest++ += a * (*x++)
						;
					})
}

template<typename T>
int arr_snlprintf(char* const s, const size_t size, size_t* lengthSoFar, const char* const preffix1,
	const char* const preffix2, const char* const singleElementFormat, const char* const separator,
	const char* const suffix1, const char* const suffix2, const T* const arr, const size_t arr_size,
	const int elementsPerLine, const bool printWithIndices, const bool indexPrintsBeforeElement)
{
	const size_t l0 = *lengthSoFar;
	int n;
	if(preffix1)
	{
		n = snlprintf(s, size, lengthSoFar, preffix1);
		if(n < 0)
			return n;
	}

	if(preffix2)
	{
		n = snlprintf(s, size, lengthSoFar, preffix2);
		if(n < 0)
			return n;
	}

	for(size_t i = 0; i < arr_size; i++)
	{
		if(i && separator)
		{
			n = snlprintf(s, size, lengthSoFar, separator);
			if(n < 0)
				return n;
		}

		if(i && elementsPerLine > 0)
		{
			if(i % elementsPerLine == 0)
			{
				n = snlprintf(s, size, lengthSoFar, "\n");
				if(n < 0)
					return n;
			}
		}

		if(!printWithIndices)
		{
			n = snlprintf(s, size, lengthSoFar, singleElementFormat, arr[i]);
			if(n < 0)
				return n;
		}
		else if(indexPrintsBeforeElement)
		{
			n = snlprintf(s, size, lengthSoFar, singleElementFormat, i, arr[i]);
			if(n < 0)
				return n;
		}
		else
		{
			n = snlprintf(s, size, lengthSoFar, singleElementFormat, arr[i], i);
			if(n < 0)
				return n;
		}
	}

	if(suffix1)
	{
		n = snlprintf(s, size, lengthSoFar, suffix1);
		if(n < 0)
			return n;
	}
	if(suffix2)
	{
		n = snlprintf(s, size, lengthSoFar, suffix2);
		if(n < 0)
			return n;
	}

	return *lengthSoFar - l0;
}

template<typename T>
int arr_printf(const char* const preffix1, const char* const preffix2, const char* const singleElementFormat,
	const char* const separator, const char* const suffix1, const char* const suffix2, const T* const arr,
	const size_t arr_size, const int elementsPerLine, const bool printWithIndices, const bool indexPrintsBeforeElement)
{
	const size_t SIZE = 4096;
	size_t len = 0;
	char S[SIZE];
	S[0] = '\0';
	int n = arr_snlprintf(S, SIZE, &len, preffix1, preffix2, singleElementFormat, separator, suffix1, suffix2, arr,
							arr_size, elementsPerLine, printWithIndices, indexPrintsBeforeElement);
	if(n < 0)
		return n;

	if(len + 1 <= SIZE)
		return fprintf(stderr, S);

	char* s = new char[n + 1];
	n = arr_snlprintf(s, n, &len, preffix1, preffix2, singleElementFormat, separator, suffix1, suffix2, arr, arr_size,
						elementsPerLine, printWithIndices, indexPrintsBeforeElement);

	if(n > 0)
		n = fprintf(stderr, s);

	delete s;
	return n;
}

template<typename T>
int arr_printf_pylist(const char* const arr_name, const char* const singleElementFormat, const T* const arr,
	const size_t arr_size)
{
	const char* const preffix2 = (arr_name) ? " = [" : "[";
	arr_printf(arr_name, preffix2, singleElementFormat, ", ", "]\n", NULL, arr, arr_size, 10, false, false);
}

template<typename T>
int arr_printf_pytuple(const char* const arr_name, const char* const singleElementFormat, const T* const arr,
	const size_t arr_size)
{
	const char* const preffix2 = (arr_name) ? " = (" : "(";
	arr_printf(arr_name, preffix2, singleElementFormat, ", ", ")\n", NULL, arr, arr_size, 10, false, false);
}

template<typename T>
int arr_printf_pydict(const char* const arr_name, const char* const singleElementFormat, const T* const arr,
	const size_t arr_size)
{
	const char* const preffix2 = (arr_name) ? " = {" : "{";
	arr_printf(arr_name, preffix2, singleElementFormat, ", ", "}\n", NULL, arr, arr_size, 10, true, true);
}

template<typename T>
int arr_printf_csv(const char* const preffix, const char* const singleElementFormat, const T* const arr,
	const size_t arr_size, const char* const separator = ", ", int elemtnsPerLine = 10)
{
	arr_printf(preffix, NULL, singleElementFormat, separator, NULL, "\n", arr, arr_size, elemtnsPerLine, false, false);
}

/*
 * Segmented array routines
 *
 * A segmented array consists of alternating segments of
 * data and undefined values.
 */
//template<typename T>
//void segarrPackByDefaultValue(T* packedarr, size_t* seglens, size_t* packedcount, size_t* segcount, size_t n, T* a,
//								T defaultValue = (T) 0)
//{
//	if (n == 0)
//	{
//		*packedcount = seglens[1] = seglens[0] = 0;
//		*segcount = 2;
//		return;
//	}
//
//	T* data2 = packedarr; // Write position in the data array.
//	size_t* segs2 = seglens; // Wirte position in the segments array.
//	size_t n1, nz; // Sizes of the data segment and default segment being packed.
//	const size_t defaultsToIgnore = sizeof(size_t) / sizeof(T);
//	const T* a_n = a + n;
//	T *a1 = a, *az = a, *a2 = a;
//	while (a2 < a_n)
//	{
//		//Begin data segment at a2.
//		for (az = a2; *az != defaultValue && az < a_n; az++)
//			;
//
//		//Begin a default segment at az.
//		for (a2 = az; *a2 == defaultValue && az < a_n; a2++)
//			;
//		//A new data segment begins at a2.
//		nz = a2 - az; // default segment length
//		if (nz > maxDefaultsToIgnore || a2 >= a_n)
//		{
//			*segs2++ = n1 = az - a1; // save the data segment length
//			*segs2++ = nz; // save the default segment length
//			for (size_t i = 0; i < n1; i++)
//				*data2++ = *a1++; // save the data from the data segment
//			a1 = a2;
//		}
//	}
//	if (nz) *segs2++ = 0;
//
//	*packedcount = data2 - packedarr;
//	*segcount = segs2 - seglens;
//}
//
//template<typename T>
//size_t segarrBufferSize(size_t packedcount, size_t segcount)
//{
//	size_t sbufsize = sizeof(size_t) * (segcount + 2);
//	size_t dbufsize = sizeof(T) * packedcount;
//	return sbufsize + dbufsize;
//}
//
//template<typename T>
//size_t segarrPackBuffer(unsigned char* buf, size_t packedcount, T* packedarr, size_t segcount, size_t* seglens)
//{
//	size_t* sbuf = (size_t*) buf;
//	size_t* sbufsize = sizeof(size_t) * segcount + 2;
//	size_t dbufsize = sizeof(T) * packedcount;
//
//	*sbuf++ = segcount;
//	*sbuf++ = seglens;
//	memcpy(sbuf, seglens, sbufsize);
//	memcpy(buf + sbufsize, packedarr, dbufsize);
//
//	return sbufsize + dbufsize;
//}
//
//template<typename T>
//void segarrUnpackBuffer(T* packedarr, size_t* packedcount, size_t* seglens, size_t* segcount, unsigned char* buf)
//{
//	size_t* sbuf = (size_t*) buf;
//	size_t* sbufsize = *sbuf++;
//	size_t* dbufsize = *sbuf++;
//	*segcount = sbufsize / sizeof(size_t) - 2;
//	*packedcount = dbufsize / sizeof(T);
//
//	memcpy(sbuf, seglens, sbufsize - 2 * sizeof(size_t));
//	memcpy(buf + sbufsize, packedarr, dbufsize);
//
//	return sbufsize + dbufsize;
//}
//
///*
// * defined in cpp
// */
////size_t segarrcounts(size_t* seglens, size_t* segcount)
//template<typename T>
//size_t segarrunpack(T* dest, T* packedarr, size_t* seglens, T defaultValue = (T) 0)
//{
//	size_t* segs = seglens;
//	if (segs[0] == 0 && segs[1] == 0) return 0;
//	T* src = packedarr;
//
//	T* dest2 = dest;
//	size_t m = *segs++;
//	while (1)
//	{
//		for (size_t j = 0; j < m; j++)
//			*dest2++ = *src++;
//		if (!(m = *segs++)) break;
//		for (size_t j = 0; j < m; j++)
//			*dest2++ = defaultValue;
//		if (!(m = *segs++)) break;
//	}
//
//	return dest2 - dest;
//}
/*
 * MPI-related
 */
#ifdef USE_MPI

template<typename T>
inline MPI_Datatype getMPI_Datatype()
{
// The following expects that the symbols corresponding to the MPI types
// being defined by the preprocessor as in OpenMPI.
#ifdef MPI_BYTE
#if __cplusplus >= 201703L /* C++17 */
	if (areOfSameType((T) 0, (byte) 0)) return MPI_BYTE;
#endif
#endif

#ifdef MPI_CHAR
	if(areOfSameType((T) 0, (char) 0))
		return MPI_CHAR;
#endif

#ifdef MPI_DOUBLE
	if(areOfSameType((T) 0, (double) 0))
		return MPI_DOUBLE;
#endif

#ifdef MPI_FLOAT
	if(areOfSameType((T) 0, (float) 0))
		return MPI_FLOAT;
#endif

#ifdef MPI_INT
	if(areOfSameType((T) 0, (int) 0))
		return MPI_INT;
#endif

#ifdef MPI_LONG
	if(areOfSameType((T) 0, (long) 0))
		return MPI_LONG;
#endif

#ifdef MPI_LONG_DOUBLE
	if(areOfSameType((T) 0, (long double) 0))
		return MPI_LONG_DOUBLE;
#endif

#ifdef MPI_LONG_INT
	if(areOfSameType((T) 0, (long int) 0))
		return MPI_LONG_INT;
#endif

#ifdef MPI_LONG_LONG
	if(areOfSameType((T) 0, (long long) 0))
		return MPI_LONG_LONG;
#endif

#ifdef MPI_LONG_LONG_INT
	if(areOfSameType((T) 0, (long long int) 0))
		return MPI_LONG_LONG_INT;
#endif

#ifdef MPI_SHORT
	if(areOfSameType((T) 0, (short) 0))
		return MPI_SHORT;
#endif

#ifdef MPI_SHORT_INT
	if(areOfSameType((T) 0, (short int) 0))
		return MPI_SHORT_INT;
#endif

#ifdef MPI_SIGNED_CHAR
	if(areOfSameType((T) 0, (signed char) 0))
		return MPI_SIGNED_CHAR;
#endif

#ifdef MPI_UNSIGNED
	if(areOfSameType((T) 0, (unsigned) 0))
		return MPI_UNSIGNED;
#endif

#ifdef MPI_UNSIGNED_CHAR
	if(areOfSameType((T) 0, (unsigned char) 0))
		return MPI_UNSIGNED_CHAR;
#endif

#ifdef MPI_UNSIGNED_LONG
	if(areOfSameType((T) 0, (unsigned long) 0))
		return MPI_UNSIGNED_LONG;
#endif

#ifdef MPI_UNSIGNED_LONG_LONG
	if(areOfSameType((T) 0, (signed long long) 0))
		return MPI_UNSIGNED_LONG_LONG;
#endif

#ifdef MPI_UNSIGNED_SHORT
	if(areOfSameType((T) 0, (unsigned short) 0))
		return MPI_UNSIGNED_SHORT;
#endif

#ifdef MPI_WCHAR
	if(areOfSameType((T) 0, (wchar_t) 0))
		return MPI_WCHAR;
#endif

#ifdef __cplusplus
#ifdef MPI_CXX_BOOL
	if(areOfSameType((T) 0, (bool) 0))
		return MPI_CXX_BOOL;
#endif
//#ifdef MPI_CXX_DOUBLE_COMPLEX
//	if (aresameargtypes((T) 0, (complex<double>) 0)) return MPI_CXX_DOUBLE_COMPLEX;
//#endif
//#ifdef MPI_CXX_FLOAT_COMPLEX
//	if (aresameargtypes((T) 0, (complex<float>) 0)) return MPI_CXX_FLOAT_COMPLEX;
//#endif
//#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
//	if (aresameargtypes((T) 0, (complex<long double>) 0)) return MPI_CXX_LONG_DOUBLE_COMPLEX;
//#endif
#else /* _cplusplus */
#ifdef MPI_C_BOOL
	if (areOfSameType((T) 0, (bool) 0)) return MPI_C_BOOL;
#endif
//#ifdef MPI_C_COMPLEX
//	if (aresameargtypes((T) 0, (complex) 0)) return MPI_C_COMPLEX;
//#endif
//#ifdef MPI_C_DOUBLE_COMPLEX
//	if (aresameargtypes((T) 0, (double complex) 0)) return MPI_C_DOUBLE_COMPLEX;
//#endif
//#ifdef MPI_C_FLOAT_COMPLEX
//	if (aresameargtypes((T) 0, (float complex) 0)) return MPI_C_FLOAT_COMPLEX;
//#endif
//#ifdef MPI_C_LONG_DOUBLE_COMPLEX
//	if (aresameargtypes((T) 0, (long double complex) 0)) return MPI_C_LONG_DOUBLE_COMPLEX;
//#endif
#endif /* __cplusplus */

//// The following C/C++ types are redundant:
////signed
////signed int
////signed long
////signed long int
////signed long long int
////signed short
////signed short int
////unsigned int
////unsigned long int
////unsigned long long
////unsigned long long int
////unsigned short int
}
template<typename T>
inline MPI_Datatype getMPI_Datatype(T a)
{
	return getMPI_Datatype<T>();
}

template<typename T>
inline MPI_Datatype getMPI_Datatype(T* a)
{
	return getMPI_Datatype<T>();
}
#endif /* USE_MPI */
/*
 template<typename T>
 void addDipoleField(T* fieldX, T* fieldY, T* fieldZ,
 const long double x, const long double y, const long double z,
 const double momentX, const double momentY, const double momentZ,
 const long double centerX, const long double centerY, const long double centerZ,
 const long double scaleFactor = 1.0)
 {
 const long double ONE_OVER_4PI = 1 / (16 * atanl(1));
 long double r, rsquared, rcube, rx, ry, rz, coeff, three_m_r_over_rsquared;
 rx = x - centerX;
 ry = y - centerY;
 rz = z - centerZ;
 rsquared = square(rx) + square(ry) + square(rz);
 rcube = cube(sqrt(rsquared));
 if(rcube)
 {
 three_m_r_over_rsquared = 3 * (momentX * rx + momentY * ry + momentZ * rz) / rsquared;
 coeff = scaleFactor * ONE_OVER_4PI / rcube;
 *fieldX += coeff * (three_m_r_over_rsquared * rx - momentX);
 *fieldY += coeff * (three_m_r_over_rsquared * ry - momentY);
 *fieldZ += coeff * (three_m_r_over_rsquared * rz - momentZ);
 }
 else
 {
 //		*fieldX = *fieldY = *fieldZ = 0;
 }
 }

 template<typename T>
 void addDipoleVectorPotential(T* vectorPotentialX, T* vectorPotentialY, T* vectorPotentialZ,
 const long double x, const long double y, const long double z,
 const double momentX, const double momentY, const double momentZ,
 const long double centerX, const long double centerY, const long double centerZ,
 const long double permeability = 1.0
 )
 {
 const long double FOUR_PI = 16 * atanl(1);
 long double r, rx, ry, rz, rsquared, coeff;
 rx = x - centerX;
 ry = y - centerY;
 rz = z - centerZ;
 rsquared = square(rx) + square(ry) + square(rz);
 if(rsquared)
 {
 coeff = permeability / FOUR_PI / rsquared
 *vectorPotentialX += coeff * (momentY * rz - momentZ * ry);
 *vectorPotentialY += coeff * (momentZ * rx - momentX * rz);
 *vectorPotentialZ += coeff * (momentX * ry - momentY * rx);
 }
 else
 {
 //		*vectorPotentialX = *vectorPotentialY = *vectorPotentialZ = 0;
 }
 }

 long double dipoleScalarPotential(const long double x, const long double y, const long double z,
 const double momentX, const double momentY, const double momentZ,
 const long double centerX, const long double centerY, const long double centerZ,
 const long double permitivity = 1.0)
 {
 const long double ONE_OVER_4PI = 1 / (16 * atanl(1));
 long double rx, ry, rz, rsquared;
 rx = x - centerX;
 ry = y - centerY;
 rz = z - centerZ;
 rsquared = square(rx) + square(ry) + square(rz);
 if(rsquared)
 {
 return ONE_OVER_4PI / permitivity / rsquared * (momentX * rx + momentY * ry + momentZ * rz);
 }
 else
 {
 return 0;
 }
 }

 template<typename T>
 void addConstFieldVectorPotential(T* vectorPotentialX, T* vectorPotentialY, T* vectorPotentialZ,
 const long double x, const long double y, const long double z,
 const double momentX, const double momentY, const double momentZ
 )
 {
 *vectorPotentialX += momentY * z;
 *vectorPotentialY += momentZ * x;
 *vectorPotentialZ += momentX * y;
 }
 */
#endif /* MY_ENZO_UTIILS_TPP */
