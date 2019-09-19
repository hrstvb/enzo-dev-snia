//#include "global_data.h"
#include <stdlib.h>
#include <stdio.h>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "Grid.h"

//int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[], float dy[],
//float dz[], FLOAT** CellLeftEdges)
//{
//	return SUCCESS;
//}

int ClearOuterVelocities(float *u, float *v, float *w, int in, int jn, int kn, int rank, float dx[], float dy[],
float dz[], FLOAT** CellLeftEdges)
{
	throw new EnzoFatalException("TODO: The Zeus solver needs to pass total energy field to ClearOuterVelocities.");
}

int ClearOuterVelocities(float *u, float *v, float *w, float *totE, int in, int jn, int kn, int rank, float dx[],
float dy[], float dz[], FLOAT** CellLeftEdges, int level, TopGridData *MetaData, int gridID)
{
//	const double TINY_V = tiny_number;
	FLOAT r_x, r_y, r_z;
	size_t index;
	bool doClearN, isOutsideR;
	double r, vdotr, rx, ry, rz, v_rx, v_ry, v_rz, v_tx, v_ty, v_tz;
	double oldVx, oldVy, oldVz, newVx, newVy, newVz;
	double oldKE, newKE;
	char negEFileName[FILENAME_MAX];
	negEFileName[0] = '\0';
	FILE *negEFlie = NULL;
	int negE1 = 0, negE2 = 0;
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
					newKE = u[index] = v[index] = w[index] = 0;
				}
				else
				{
					// Project the velocity onto the radius.
					vdotr = r_x * u[index] + r_y * v[index] + r_z * w[index];

					// Determine if the radial component should be cleared
					// based on the parameters and depending on whether it is
					// pointing toward the center (inward) or away (outward).
					doClearN = (OuterVelocitiesClearInward && vdotr < 0) || (OuterVelocitiesClearOutward && vdotr > 0);

					if(OuterVelocitiesClearTangential && doClearN)
					{
						// If clear both the tangential and the radial components,
						// set the velocity to zero and move to the next zone.
						newKE = u[index] = v[index] = w[index] = 0;
					}
					else if(!OuterVelocitiesClearTangential && !doClearN)
					{
						// If clear neither tangential nor the radial component,
						// there is nothing to do for this zone, so move to the next one.
						newKE = oldKE;
					}
					else
					{
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
							newVx = oldVx - v_rx;
							newVy = oldVy - v_ry;
							newVz = oldVz - v_rz;
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
					}
				}

				float dE = newKE - oldKE;
				float oldE = totE[index];
				float newE = oldE + dE;
				int fixType = 0;
				if(oldE <= 0)
				{
					TRACEF("Total energy <=0 before correction   %e    at    %e  %e  %e  ( %lld  %lld  %lld )",
							oldE, r_x, r_y, r_z, i, j, k);
					newE = (oldKE > 0) ? oldKE : tiny_pressure / (1 - Gamma);
					fixType = 1;
					negE1++;
				}
				else if(newE <= 0)
				{
					TRACEF("Corrected total energy <= 0   %e    at    %e  %e  %e  ( %lld  %lld  %lld )", newE,
							r_x, r_y, r_z, i, j, k);
					newE = tiny_pressure / (1 - Gamma);
					fixType = 2;
					negE2++;
				}
				totE[index] = newE;
				if(fixType)
				{
					if(negEFlie == NULL)
					{
						char cmd[FILENAME_MAX + 10];
						char filedir[FILENAME_MAX + 1];
						cmd[0] = filedir[0] = '\0';
						if(MetaData->GlobalDir)
							snprintf(filedir, FILENAME_MAX, "%s/", MetaData->GlobalDir);
						else if(MetaData->LocalDir)
							snprintf(filedir, FILENAME_MAX, "%s/", MetaData->LocalDir);

						snprintf(filedir, FILENAME_MAX, "%snege/ne%06lld", filedir, MetaData->CycleNumber);
						snprintf(cmd, FILENAME_MAX + 9, "mkdir -p %s", filedir);
						system(cmd);
						snprintf(negEFileName, FILENAME_MAX, "%s/ne%06lld_%lld_%02lld.py", filedir,
									MetaData->CycleNumber, level, gridID);
						negEFlie = fopen(negEFileName, "w+");
						fprintf(negEFlie, "# fixType, i, j, k\n"
								"# fixType = 1 -- total energy <= 0 before subtracting kinetic energy;\n"
								"# fixType = 2 -- total energy <= 0 after subtracting kinetic energy;\n"
								"#\n"
								"type('', (), dict(\n"
								"cycle = %f,\n"
								"time = %f,\n"
								"cols = 'level gridID fixType i j k'.split(),\n"
								"level = %lld,\n"
								"gridID = %lld,\n"
								"gridDims = ( %lld , %lld , %lld ),\n"
								"gridLeftEdge = ( %e , %e , %e ),\n"
								"gridRightEdge = ( %e , %e , %e ),\n"
								"nGhostZones = %lld,\n"
								"data = np.array([\n",
								MetaData->CycleNumber, MetaData->Time, MetaData->DataDumpDir, level, gridID, in, jn, kn,
								CellLeftEdges[0][0], CellLeftEdges[1][0], CellLeftEdges[2][0], CellLeftEdges[0][in],
								CellLeftEdges[1][jn], CellLeftEdges[2][kn], NumberOfGhostZones);
					}
					fprintf(negEFlie, "( %lld , %lld , %lld , %lld , %lld , %lld ),", fixType, i, j, k);
				}
			}
		}
	}

	if(negE1 + negE2)
	{
		TRACEF(" Level %lld, grid %lld: negative energy in %lld zones (%lld before, %lld after correction)", level,
				gridID, negE1, negE1, negE2);
	}
	if(negEFlie && negEFlie != stderr && negEFlie != stdout)
	{
		fprintf(negEFlie, "])))\n");
		fclose(negEFlie);
		TRACEF("'%s' written.\n", negEFileName);
	}

	return SUCCESS;
}

int grid::ClearOuterVelocities(int level, TopGridData *MetaData)
{
	if(MyProcessorNumber != ProcessorNumber)
		return SUCCESS;

	float *vxField, *vyField, *vzField, *totEField;
	MHD_SNIA_GetFields(NULL, &totEField, NULL, &vxField, &vyField, &vzField, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	NULL);

	return ::ClearOuterVelocities(vxField, vyField, vzField, totEField, GridDimension[0], GridDimension[1],
									GridDimension[2], GridRank, CellWidth[0], CellWidth[1], CellWidth[2], CellLeftEdge,
									level, MetaData, this->ID);
}
