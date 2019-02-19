/***********************************************************************
 /
 /  GRID CLASS (Compute Del^2 on a field)
 /
 /  written by:  Boyan Hristov
 /  date:        August, 2014
 /
 /  PURPOSE:  Computes the Laplacian of 'sourceField' and returns
 /		the result in the 'resultField'.
 /
 /	'method' determines the numerical alhorithm:
 /    	0 = finite difference using a 3 (1D), 5 (2D) or 7 (3D) point stencil;
 /		1 = finite difference using a 27 point stencil in 3D;
 /		2 = finite difference using a 125 point stencil in 3D.
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/
#include <math.h>
#include <stdio.h>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

int grid::ComputeLaplacian(float* resultField, float* sourceField, int method)
{
	int offset, di, dj, dk;
	float *result, *source;
	FLOAT dxdx = POW(CellWidth[0][0], 2);
	//TODO: dx *= LengthUnits?

	//if ( !isCartesian )
	//	ENZO_FAIL("grid::ComputeLaplacian: grid is not cartesian.");

	switch(method)
	{
	case 0:
		switch(GridRank)
		{
		case 1:

			offset = GridStartIndex[0];
			source = sourceField + offset;
			result = resultField + offset;

			for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
			{
				*result = (*(source + 1) + *(source - 1) - 2 * (*source)) / dxdx;

				source++;
				result++;

			} // for i

			break;

		case 2:
			di = ELT(1, 0, 0); // - ELT( 0, 0, 0 );
			dj = ELT(0, 1, 0); // - ELT( 0, 0, 0 );

			for(int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
			{
				offset = ELT(GridStartIndex[0], j, 0);
				source = sourceField + offset;
				result = resultField + offset;

				for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
				{
					*result = (*(source + di) + *(source - di) + *(source + dj) + *(source - dj) - 4 * (*source))
							/ dxdx;

					source++;
					result++;

				} // for i
			} // for j

			break;

		case 3:
			di = ELT(1, 0, 0); // - ELT( 0, 0, 0 );
			dj = ELT(0, 1, 0); // - ELT( 0, 0, 0 );
			dk = ELT(0, 0, 1); // - ELT( 0, 0, 0 );

			for(int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
			{
				for(int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
				{
					offset = ELT(GridStartIndex[0], j, k);
					source = sourceField + offset;
					result = resultField + offset;

					for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
					{
						*result = (*(source + di) + *(source - di) + *(source + dj) + *(source - dj) + *(source + dk)
								+ *(source - dk) - 6 * (*source)) / dxdx;

						source++;
						result++;

					} // for i
				} // for j
			} // for k

			break;
		} // end switch GridRank

		break;

		//end switch method 0

	case 1:
		switch(GridRank)
		{
		case 1:

			offset = GridStartIndex[0];
			source = sourceField + offset;
			result = resultField + offset;

			for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
			{
				*result = (*(source + 1) + *(source - 1) - 2 * (*source)) / dxdx;

				source++;
				result++;

			} // for i

			break;

		case 3:
			di = ELT(1, 0, 0); // - ELT( 0, 0, 0 );
			dj = ELT(0, 1, 0); // - ELT( 0, 0, 0 );
			dk = ELT(0, 0, 1); // - ELT( 0, 0, 0 );

			for(int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
			{
				for(int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
				{
					offset = ELT(GridStartIndex[0], j, k);
					source = sourceField + offset;
					result = resultField + offset;

					for(int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
					{
						*result = ( *(source + di + dj + dk)
							  + *(source + di + dj +  0)
							  + *(source + di + dj - dk)
							  + *(source + di +  0 + dk)
							  + *(source + di +  0 +  0)
							  + *(source + di +  0 - dk)
							  + *(source + di - dj + dk)
							  + *(source + di - dj +  0)
							  + *(source + di - dj - dk)

							  + *(source +  0 + dj + dk)
							  + *(source +  0 + dj +  0)
							  + *(source +  0 + dj - dk)
							  + *(source +  0 +  0 + dk)
							  //+ *(source +  0 +  0 +  0)
							  + *(source +  0 +  0 - dk)
							  + *(source +  0 - dj + dk)
							  + *(source +  0 - dj +  0)
							  + *(source +  0 - dj - dk)

							  + *(source - di + dj + dk)
							  + *(source - di + dj +  0)
							  + *(source - di + dj - dk)
							  + *(source - di +  0 + dk)
							  + *(source - di +  0 +  0)
							  + *(source - di +  0 - dk)
							  + *(source - di - dj + dk)
							  + *(source - di - dj +  0)
							  + *(source - di - dj - dk)

							  - 26 * (*source)

							  ) / 9 / dxdx;

						source++;
						result++;

					} // for i
				} // for j
			} // for k

			break;
		} // end switch GridRank

		break;
	case 2:
		switch(GridRank)
		{
		case 3:
//			-88 0   0   0
//			4   0   1   1
//			4   0   2   2
//			1   0   3   3
//			-2   1   0   4
//			1   1   1   5
//			1   1   2   6
//			-1   2   0   8
//			-1   2   1   9
//			1   3   0   12
//			/24

			// ones                  0  1  2  3   0  1  2  -   0   1  -  -  0
			// twos                  0  1  2  3   1  1  1  -   2   2  -  -  3
			// ones^2 + twos^2       0  1  2  3   4  5  6  -   8   9  -  - 12
			const int cR[13] = {-88, 4, 4, 1, -2, 1, 1, 0, -1, -1, 0, 0, 1};
			const double cD = 24;

			const size_t i1 = GridStartIndex[0];
			const size_t i2 = GridEndIndex[0];
			const size_t j1 = GridStartIndex[1];
			const size_t j2 = GridEndIndex[1];
			const size_t k1 = GridStartIndex[2];
			const size_t k2 = GridEndIndex[2];
			const size_t n = i2 - i1 + 1;

			int cIJK[125];
			int* c = cIJK;
			int ci = 0, RR;

			for(int KK = -2; KK <= 2; KK++)
				for(int JJ = -2; JJ <= 2; JJ++)
					for(int II = -2; II <= 2; II++)
					{
						RR = II * II + JJ * JJ + KK * KK;
						cIJK[ci++] = cR[RR];
					}

			arr_set(resultField, GetGridSize(), 0);

			for(size_t k = k1; k <= k2; k++)
				for(size_t j = j1; j <= j2; j++)
				{
					result = resultField + ELT(i1, j, k);
					c = cIJK;
					for(int KK = -2; KK <= 2; KK++)
						for(int JJ = -2; JJ <= 2; JJ++)
							for(int II = -2; II <= 2; II++)
							{
								source = sourceField + ELT(i1 + II, j + JJ, k + KK);
								arr_axpy(result, source, n, *c++);
							}
				}

			break;
		}
		break;
	} //end switch method 2

	return SUCCESS;
}

