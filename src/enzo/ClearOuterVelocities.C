//#include "global_data.h"
#include "myenzoutils.h"
#include "ErrorExceptions.h"
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

int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[],
float dy[], float dz[], FLOAT** CellLeftEdges)
{
	throw new EnzoFatalException("TODO: The Zeus solver needs to pass total energy field to ClearOuterVelocities.");
}

int ClearOuterVelocities(float *u, float *v, float *w, float *totE, int in, int jn, int kn, int rank, float dx[],
float dy[], float dz[], FLOAT** CellLeftEdges)
{
//	const double TINY_V = tiny_number;
	FLOAT r_x, r_y, r_z;
	size_t index;
	bool doClearN, isOutsideR;
	double r, vdotr, rx, ry, rz, v_rx, v_ry, v_rz, v_tx, v_ty, v_tz;
	double oldVx, oldVy, oldVz, newVx, newVy, newVz;
	double oldKE, newKE;
	int flagsIOT = 4 * OuterVelocitiesClearInward + 2 * OuterVelocitiesClearOutward + OuterVelocitiesClearTangential;

	if(OuterVelocitiesSphereRadius < 0 || flagsIOT == 0)
		return SUCCESS;

	v_rx = v_ry = v_rz = 0;
	// Clear velocity beyond a certain radius
	for(int k = 0; k < kn; k++)
	{
		for(int j = 0; j < jn; j++)
		{
			for(int i = 0; i < in; i++)
			{
				r_x = CellLeftEdges[0][i] - SphericalGravityCenter[0];
				r_y = CellLeftEdges[1][j] - SphericalGravityCenter[1];
				r_z = CellLeftEdges[2][k] - SphericalGravityCenter[2];

				if(HydroMethod == Zeus_Hydro)
				{
					isOutsideR = false;
					isOutsideR |= lenl(r_x - dx[i] / 2, r_y, r_z) > OuterVelocitiesSphereRadius;
					isOutsideR |= lenl(r_x, r_y - dy[i] / 2, r_z) > OuterVelocitiesSphereRadius;
					isOutsideR |= lenl(r_x, r_y, r_z - dz[i] / 2) > OuterVelocitiesSphereRadius;
				}
				else
				{
					r = lenl(r_x, r_y, r_z);
					isOutsideR = r > OuterVelocitiesSphereRadius;
				}

				// Do not change velocities inside the radius.
				if(!isOutsideR)
					continue;

				index = i + in * (j + jn * k);
				oldVx = u[index];
				oldVy = v[index];
				oldVz = w[index];
				oldKE = (square(oldVx) + square(oldVy) + square(oldVz)) / 2;

				// If clear the normal component regardless of direction
				// (i.e both inward and outward) and clear the tangential
				// then set velocity to zero and move to the next cell.
				if(flagsIOT == 7)
				{
					u[index] = v[index] = w[index] = 0;
					totE[index] -= oldKE;
					if(totE[index] < 0)
						TRACEF("Negative energy    %e    at    %e  %e  %e  ( %lld  %lld  %lld )", totE[index], r_x, r_y,
								r_z, i, j, k);
					continue;
				}

				// Project the velocity onto the radius.
				vdotr = r_x * u[index] + r_y * v[index] + r_z * w[index];

				// Determine if the radial component should be cleared
				// based on the parameters and depending on whether it is
				// pointing toward the center (inward) or away (outward).
				doClearN = (OuterVelocitiesClearInward && vdotr < 0) || (OuterVelocitiesClearOutward && vdotr > 0);

				if(OuterVelocitiesClearTangential)
				{
					// If clear both the tangential and the radial components,
					// set the velocity to zero and move to the next zone.
					if(doClearN)
					{
						u[index] = v[index] = w[index] = 0;
						totE[index] -= oldKE;
						if(totE[index] < 0)
							TRACEF("Negative energy    %e    at    %e  %e  %e  ( %lld  %lld  %lld )", totE[index], r_x,
									r_y, r_z, i, j, k);
						continue;
					}
				}
				else
				{
					// If clear neither tangential nor the radial component,
					// there is nothing to do for this zone, so move to the next one.
					if(!doClearN)
						continue;
				}

				// We need to clear either the radial or the tangential
				// component of the velocity and keep the other.

				if(HydroMethod == Zeus_Hydro)
				{
					r = lenl(r_x, r_y, r_z);
				}

				// Get radial velocity:
				r *= r;
				v_rx = vdotr * r_x / r;
				v_ry = vdotr * r_y / r;
				v_rz = vdotr * r_z / r;

				if(doClearN)
				{
//					TRACEF(" %e %e %e    %e %e %e    %e %e %e    %e %e", u[index], v[index], w[index], v_rx, v_ry, v_rz, r_x, r_y, r_z, vdotr, r);
					newVx = oldVx - v_rx;
					newVy = oldVy - v_ry;
					newVz = oldVz - v_rz;
//					TRACEF(" %e %e %e    %e %e %e    %e %e %e    %e %e", u[index], v[index], w[index], v_rx, v_ry, v_rz, r_x, r_y, r_z, vdotr, r);
				}
				else
				{
					newVx = v_rx;
					newVy = v_ry;
					newVz = v_rz;
				}

				newKE = (square(newVx) + square(newVy) + square(newVz)) / 2;

				u[index] = newVx;
				v[index] = newVy;
				w[index] = newVz;
				totE[index] += newKE - oldKE;
				if(totE[index] < 0)
					TRACEF("Negative energy    %e    at    %e  %e  %e  ( %lld  %lld  %lld )", totE[index], r_x, r_y,
							r_z, i, j, k);
			}
		}
	}
	return SUCCESS;
}

int grid::ClearOuterVelocities()
{
	if(MyProcessorNumber != ProcessorNumber)
		return SUCCESS;

	float *vxField, *vyField, *vzField, *totEField;
	MHD_SNIA_GetFields(NULL, &totEField, NULL, &vxField, &vyField, &vzField, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	NULL);

	return ::ClearOuterVelocities(vxField, vyField, vzField, totEField, GridDimension[0], GridDimension[1],
									GridDimension[2], GridRank, CellWidth[0], CellWidth[1], CellWidth[2], CellLeftEdge);
}
