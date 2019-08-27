//#include "global_data.h"
#include "myenzoutils.h"
#include "Grid.h"

//int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[], float dy[],
//float dz[], FLOAT** CellLeftEdges)
//{
//	if(OuterVelocitiesSphereRadius < 0 || (!OuterVelocitiesClearInward && !OuterVelocitiesClearOutward))
//		return SUCCESS;
//
//	const double TINY_V = tiny_number;
//	FLOAT x, y, z;
//	size_t index;
//	bool doClean1, doClean2;
//	double rv;
//
//	// Clear velocity beyond a certain radius
//	//for(int k = 0; k < kn && ZEUS_IncludeDivergenceTerm; k++)
//	for(int k = 0; k < kn; k++)
//	{
//		for(int j = 0; j < jn; j++)
//		{
//			for(int i = 0; i < in; i++)
//			{
//				x = CellLeftEdges[0][i] - SphericalGravityCenter[0];
//				y = CellLeftEdges[1][j] - SphericalGravityCenter[1];
//				z = CellLeftEdges[2][k] - SphericalGravityCenter[2];
//				index = IDX(i, j, k);
//
//				double rv = x * u[index] + y * v[index] + z * w[index];
//				doClean1 = false;
//				doClean1 |= OuterVelocitiesClearInward && rv < -TINY_V;
//				doClean1 |= OuterVelocitiesClearOutward && rv > TINY_V;
//				doClean1 |= OuterVelocitiesClearTangential && (-TINY_V <= rv) && (rv <= TINY_V);
//
//				if(HydroMethod == Zeus_Hydro)
//				{
//					doClean2 = false;
//					doClean2 |= lenl(x - dx[i] / 2, y, z) > OuterVelocitiesSphereRadius;
//					doClean2 |= lenl(x, y - dy[i] / 2, z) > OuterVelocitiesSphereRadius;
//					doClean2 |= lenl(x, y, z - dz[i] / 2) > OuterVelocitiesSphereRadius;
//				}
//				else
//				{
//					doClean2 = lenl(x, y, z) > OuterVelocitiesSphereRadius;
//				}
//
//				if(doClean1 && doClean2)
//				{
//					u[index] = v[index] = w[index] = 0;
//				}
//			}
//		}
//	}
//	return SUCCESS;
//}

int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[], float dy[],
float dz[], FLOAT** CellLeftEdges)
{
	if(OuterVelocitiesSphereRadius < 0
			|| (!OuterVelocitiesClearInward && !OuterVelocitiesClearOutward && !OuterVelocitiesClearTangential))
		return SUCCESS;

//	const double TINY_V = tiny_number;
	FLOAT x, y, z;
	size_t index;
	bool doClean1, doClean2;
	double r, vdotr, rx, ry, rz, v_rx, v_ry, v_rz, v_tx, v_ty, v_tz;

	// Clear velocity beyond a certain radius
	if(OuterVelocitiesClearInward && OuterVelocitiesClearOutward && OuterVelocitiesClearTangential)
	{
		for(int k = 0; k < kn; k++)
		{
			for(int j = 0; j < jn; j++)
			{
				for(int i = 0; i < in; i++)
				{
					x = CellLeftEdges[0][i] - SphericalGravityCenter[0];
					y = CellLeftEdges[1][j] - SphericalGravityCenter[1];
					z = CellLeftEdges[2][k] - SphericalGravityCenter[2];
					index = i + in * (j + jn * k);

					if(HydroMethod == Zeus_Hydro)
					{
						doClean2 = false;
						doClean2 |= lenl(x - dx[i] / 2, y, z) > OuterVelocitiesSphereRadius;
						doClean2 |= lenl(x, y - dy[i] / 2, z) > OuterVelocitiesSphereRadius;
						doClean2 |= lenl(x, y, z - dz[i] / 2) > OuterVelocitiesSphereRadius;
					}
					else
					{
						doClean2 = lenl(x, y, z) > OuterVelocitiesSphereRadius;
					}

					if(doClean2)
					{
						u[index] = v[index] = w[index] = 0;
					}
				}
			}
		}
		return SUCCESS;
	}

	if(OuterVelocitiesClearInward && (OuterVelocitiesClearOutward || OuterVelocitiesClearTangential))
	{
		for(int k = 0; k < kn; k++)
		{
			for(int j = 0; j < jn; j++)
			{
				for(int i = 0; i < in; i++)
				{
					x = CellLeftEdges[0][i] - SphericalGravityCenter[0];
					y = CellLeftEdges[1][j] - SphericalGravityCenter[1];
					z = CellLeftEdges[2][k] - SphericalGravityCenter[2];
					index = i + in * (j + jn * k);

					vdotr = x * u[index] + y * v[index] + z * w[index];
					doClean1 = (OuterVelocitiesClearOutward) ? (vdotr >= 0) : (vdotr <= 0);

					if(HydroMethod == Zeus_Hydro)
					{
						doClean2 = false;
						doClean2 |= lenl(x - dx[i] / 2, y, z) > OuterVelocitiesSphereRadius;
						doClean2 |= lenl(x, y - dy[i] / 2, z) > OuterVelocitiesSphereRadius;
						doClean2 |= lenl(x, y, z - dz[i] / 2) > OuterVelocitiesSphereRadius;
					}
					else
					{
						doClean2 = lenl(x, y, z) > OuterVelocitiesSphereRadius;
					}

					if(doClean2)
					{
						u[index] = v[index] = w[index] = 0;
					}
				}
			}
		}
		return SUCCESS;
	}

	v_rx = v_ry = v_rz = v_tx = v_ty = v_tz = 0;
	for(int k = 0; k < kn; k++)
	{
		for(int j = 0; j < jn; j++)
		{
			for(int i = 0; i < in; i++)
			{
				x = CellLeftEdges[0][i] - SphericalGravityCenter[0];
				y = CellLeftEdges[1][j] - SphericalGravityCenter[1];
				z = CellLeftEdges[2][k] - SphericalGravityCenter[2];
				index = i + in * (j + jn * k);

				r = lenl(x, y, z);

				if(HydroMethod == Zeus_Hydro)
				{
					doClean2 = false;
					doClean2 |= lenl(x - dx[i] / 2, y, z) > OuterVelocitiesSphereRadius;
					doClean2 |= lenl(x, y - dy[i] / 2, z) > OuterVelocitiesSphereRadius;
					doClean2 |= lenl(x, y, z - dz[i] / 2) > OuterVelocitiesSphereRadius;
				}
				else
				{
					doClean2 = r > OuterVelocitiesSphereRadius;
				}

				if(!doClean2)
					continue;

				if(r < tiny_number)
					continue;

				rx = x / r;
				ry = y / r;
				rz = z / r;

				vdotr = rx * u[index] + ry * v[index] + rz * w[index];
				v_rx = rx * vdotr;
				v_ry = ry * vdotr;
				v_rz = rz * vdotr;
				if(OuterVelocitiesClearTangential)
				{
					doClean1 = true;
				}
				else
				{
					doClean1 = false;
					v_tx = u[index] - v_rx;
					v_ty = u[index] - v_ry;
					v_tz = u[index] - v_rz;
				}
				if(OuterVelocitiesClearInward)
					doClean1 |= vdotr < 0;
				if(OuterVelocitiesClearOutward)
					doClean1 |= vdotr > 0;

				if(doClean1)
				{
					u[index] = v_rx + v_tx;
					v[index] = v_ry + v_ty;
					w[index] = v_rz + v_tz;
				}
			}
		}
	}
	return SUCCESS;
}

int grid::ClearOuterVelocities()
{
	float *vxField, *vyField, *vzField;
	MHD_SNIA_GetFields(NULL, NULL, NULL, &vxField, &vyField, &vzField, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	NULL);
	return ::ClearOuterVelocities(vxField, vyField, vzField, GridDimension[0], GridDimension[1], GridDimension[2],
									GridRank, CellWidth[0], CellWidth[1], CellWidth[2], CellLeftEdge);
}
