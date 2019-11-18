#include <cmath>
#include <float.h>

#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "mylimiters.h"
#include "global_data.h"

#define MYLIM_FORMAT ".15f"
#define MYLIM_EPSILON (DBL_EPSILON)
#define MYLIM_TOLL (MYLIM_EPSILON*16)

inline MYLIM_FLOAT limiter1n(const int n, const MYLIM_FLOAT x)
{
	MYLIM_FLOAT xn = pow(x, n);
	return 2 * (square(xn) + x) / square(xn + 1);
}

inline MYLIM_FLOAT limiter12(const int n, const MYLIM_FLOAT x)
{
	MYLIM_FLOAT xx = square(x);
	return 2 * (square(xx) + x) / square(xx + 1);
}

inline MYLIM_FLOAT limiter11(const int n, const MYLIM_FLOAT x)
{
	return 2 * x / (x + 1);
}

inline MYLIM_FLOAT limiter(const int n, MYLIM_FLOAT r)
{
	switch(n)
	{
	case 1:
	case 11:
		return limiter11(n, r);
	case 2:
	case 12:
		return limiter12(n, r);
	default:
		return limiter1n(n, r);
	}
}

MYLIM_FLOAT invLimiter(int n, MYLIM_FLOAT dP01, MYLIM_FLOAT slope, int *err)
{
	bool debug = *err;
	bool upperBound = false;
	*err = 0;
	int sign01 = sign<MYLIM_FLOAT>(dP01);
	int sign02 = sign<MYLIM_FLOAT>(slope);

	if(sign01 * sign02 <= 0)
	{
		*err = 1;
		return 0;
	}

	slope /= dP01;
	if(slope == 0)
	{
		*err = 2;
		return 0;
	}

	if(slope == 1)
	{
		*err = 3;
		return dP01;
	}

	// slope > 0 sould be guaranteed by the above code.
	//if(slope < 0)
	//{
	//	*err = -1;
	//	ENZO_VFAIL("Negative slope, %f (%e) for the inverse of MHDCTSlopeLimiter=%lld.\n", slope, slope, n + 10);
	//}

	if(slope >= 2)
	{
		*err = -2;
		if(upperBound)
		{
			TRACEF("Greater than 2 slope, %f (%e) for the inverse of MHDCTSlopeLimiter=%lld.\n", slope, slope, n + 10);
			throw new EnzoFatalException("", __FILE__, __LINE__);
//			ENZO_VTHROW("Greater than 2 slope, %f (%e) for the inverse of MHDCTSlopeLimiter=%lld.\n", slope, slope,
//						n + 10);
		}
	}

	double l, m, r, y;
	double (*limfn)(int, double);

	switch(n)
	{
	case 1:
		return dP01 * abs(slope / (2 * dP01 - slope));
	case 2:
		limfn = limiter12;
		break;
	default:
		if(n < 1 || n > 9)
		{
			TRACEF("Inverse for MHDCTSlopeLimiter=%lld not implemented.\n", n + 10);
			throw new EnzoFatalException("", __FILE__, __LINE__);
			//ENZO_VTHROW("Inverse for MHDCTSlopeLimiter=%lld not implemented.\n", n + 10);
		}

		if(slope < MyLimiterX12Y12[n][2])
		{
			limfn = limiter1n;
			*err = 4;
			break;
		}
		else if(MyLimiterX12Y12[n][3] < slope && upperBound)
		{
			limfn = limiter1n;
			*err = 5;
			break;
		}
		else
		{
			return dP01 * (2 * slope - 1);
		}
		limfn = limiter1n;
	}

	const MYLIM_FLOAT TOLL = DBL_EPSILON * 16;
	if(slope < 1)
	{
		l = 0;
		r = MyLimiterX12Y12[n][0];
	}
	else
	{
		l = MyLimiterX12Y12[n][1];
		r = 4;
		while(1)
		{
			y = (*limfn)(n, r);
//			if(debug)
//				TRACEF("-------- l,m,r= %e %e %e y,s,y-s= %e %e %e r-l= %e", l, m, r, y, slope, y-slope, r-l);
			switch(sign<double>(y - slope))
			{
			case -1:
				r *= 2;
				continue;
			case 0:
				l = m = r;
			}
			break;
		}
	}

	// Find the root between l and r.
	while(r - l > TOLL)
	{
		m = (l + r) / 2;
		y = (*limfn)(n, m);
//			if(debug)
//			TRACEF("------------------- l,m,r= %e %e %e y,s,y-s= %e %e %e r-l= %e", l, m, r, y, slope, y-slope, r-l);

		switch(sign<MYLIM_FLOAT>(y - slope))
		{
		case -1:
			l = m;
			continue;
//		case 0:
//			break;
		default:
			r = m;
			continue;
		}
		break;
	}

	return m * dP01;
}

MYLIM_FLOAT invLimiter(MYLIM_FLOAT dP01, MYLIM_FLOAT slope, int *err)
{
	int n = MHDCTSlopeLimiter - 10;
	return invLimiter(n, dP01, slope, err);
}

void limiterTest()
{
	if(MyProcessorNumber != ROOT_PROCESSOR)
		return;
	TRACEF("LIMITERS TEST BEGIN ----------------------------");
	const int N = 20;
	for(int n = 3; n < 4; n++)
	{
		for(int i = 1; i < N; i++)
		{
			int err = 0;
			MYLIM_FLOAT s1 = (MYLIM_FLOAT) 2 * (MYLIM_FLOAT) i / (MYLIM_FLOAT) N;
			MYLIM_FLOAT r = invLimiter(n, 1, s1, &err);
			MYLIM_FLOAT s2 = limiter(n, r);
			TRACEF("n,i= %lld %lld s1,x,s2= %e %e %e s1-s2,s2/s1-1= %e %e", n, i, s1, r, s2, s2 - s1, s2 / s1 - 1);
		}
	}
	TRACEF("LIMITERS TEST END ------------------------------");
}

void limiter1nInit()
{
	if(MyLimiterX12Y12[0][0])
		return;

	const MYLIM_FLOAT TOLL = DBL_EPSILON * 16;
	arr_set(MyLimiterX12Y12[0], 40, 1);

	for(int n = 3; n < 10; n++)
	{
		MYLIM_FLOAT l, m, r, y, s;
		for(int i = 0; i < 2; i++)
		{
			m = l = (i) ? 2 : 0.25;
			r = (i) ? 4 : 0.5;
			s = (i) ? (-1) : 1;

			while(r - l > TOLL)
			{
				m = (l + r) / 2;
				y = limiter1n(n, m);

				switch(sign<MYLIM_FLOAT>(s * (y - (m + 1) / 2)))
				{
				case -1:
					l = m;
					continue;
				case 0:
					break;
				default:
					r = m;
					continue;
				}
				break;
			}

			MyLimiterX12Y12[n][i] = m;
			MyLimiterX12Y12[n][i + 2] = y;
		}
	}

	if(MyProcessorNumber == ROOT_PROCESSOR)
		arr_printf_csv("SlopeLimiter X1, X2, Y1, Y2 = \n", "%" MYLIM_FORMAT, MyLimiterX12Y12[0], 40, " ", 4);

	limiterTest();
}
