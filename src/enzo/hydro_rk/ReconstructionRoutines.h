//#include "macros_and_parameters.h"
//#include <math.h>
//#include <stdio.h>

inline float Max(float a1, float a2, float a3, float a4, float a5, float a6,
	  float a7, float a8, float a9, float a10, float a11, float a12)
{
  float a = a1;

  if (a2 > a)
    a = a2;

  if (a3 > a)
    a = a3;

  if (a4 > a)
    a = a4;

  if (a5 > a)
    a = a5;

  if (a6 > a)
    a = a6;

  if (a7 > a)
    a = a7;

  if (a8 > a)
    a = a8;

  if (a9 > a)
    a = a9;

  if (a10 > a)
    a = a10;

  if (a11 > a)
    a = a11;

  if (a12 > a)
    a = a12;

  return a;
}

inline float Max(float a1, float a2, float a3, float a4, float a5, float a6,
		    float a7, float a8)
{
  float a = a1;

  if (a2 > a)
    a = a2;

  if (a3 > a)
    a = a3;

  if (a4 > a)
    a = a4;

  if (a5 > a)
    a = a5;

  if (a6 > a)
    a = a6;

  if (a7 > a)
    a = a7;

  if (a8 > a)
    a = a8;

  return a;
}


inline float Max(float a, float b, float c)
{
  if (a > b) {
    if (a > c)
      return a;
    else
      return c;
  } else {
    if (b > c)
      return b;
    else
      return c;
  }
}

inline float Min(float a, float b, float c)
{
  if (a<b) {
    if (c<a)
      return c;
    else
      return a;
  } else {
    if (c<b)
      return c;
    else
      return b;
  }
}


inline float minmod(float a, float b)
{
  /*if (a*b <= 0) {
    return 0.0;
    } */
  return 0.5*(sign(a)+sign(b))*min(fabs(a), fabs(b));
}


inline float minmod(float a, float b, float c)
{

  /*if (a*b <=0)
    return 0.0;
  if (b*c <=0)
  return 0.0;*/

  /*if (r1 != r2) {
    printf("r1!=r2: r1 = %lf, r2 = %lf, sign(a)=%d, sign(b)=%d, sign(c)=%d\n",
	   r1, r2, sign(a), sign(b), sign(c));

  }
  return r1;*/
  //return sign(a)*Min(fabs(a), fabs(b), fabs(c));
  return 0.25*(sign(a)+sign(b))*ABS((sign(a)+sign(c)))*Min(fabs(a), fabs(b), fabs(c));
}

inline float minmod(float a, float b, float c, float *minmodtwin, float a2, float b2, float c2)
{
	//return 0.25*(sign(a)+sign(b))*ABS((sign(a)+sign(c)))*Min(fabs(a), fabs(b), fabs(c));
//	*minmodtwin = c2; // central difference;
//	*minmodtwin = minmod(a2, b2, c2); // regular minmod;
//	return minmod(a, b, c);

	// Twin minmod:
	int s = sign(a);
	s = (s + sign(b)) * ABS(s + sign(c)) / 4;
	if(s == 0)
		return *minmodtwin = 0;

	a = fabs(a);
	b = fabs(b);
	c = fabs(c);

	//min(a,b,c)=?
	if(a<b)
	{
		if(a<c)
		{
			// a<b,c, min=a
			*minmodtwin = a2;
			return s * a;
		}
		//c<a<b
	}
	else //i.e. b<a
	{
		if(b<c)
		{
			//b<a,c, min=b
			*minmodtwin = b2;
			return s * b;
		}
		// c<b<a
	}

	//c<a,b, min=c
	*minmodtwin = c2;
	return s * c;
}

inline int same_sign(float a, float b, float c)
{
  /*
  if (a*b <= 0)
    return 0;
  if (b*c <= 0)
    return 0;
    return 1;*/
  return 0.25*(sign(a)+sign(b))*(sign(a)+sign(c));
}

inline int which_MIN(float a, float b, float c)
{
  if (a < b) {
    if (c < a)
      return 2;
    else
      return 0;
  } else {
    if (c < b)
      return 2;
    else
      return 1;
  }

}


