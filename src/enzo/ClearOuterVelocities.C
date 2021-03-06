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

FILE *negEFile_open(char** filename, TopGridData *MetaData, int level, grid *g, char* dumpPrefix, char* dumpSuffix,
	char* dumpPreamble1, char* dumpPreamble2)
{
	FILE *file;
	char cmd[FILENAME_MAX + 10];
	char dirname[FILENAME_MAX + 1];
	char *fname = (filename && *filename) ? *filename : new char[FILENAME_MAX + 1];
	size_t n;

	cmd[0] = dirname[0] = fname[0] = '\0';
	if(MetaData->GlobalDir)
		n = snprintf(dirname, FILENAME_MAX, "%s", MetaData->GlobalDir);
	else if(MetaData->LocalDir)
		n = snprintf(dirname, FILENAME_MAX, "%s", MetaData->LocalDir);

	if(n && dirname[n] != '/')
		n += snprintf(dirname + n, FILENAME_MAX, "/");

	n += snprintf(dirname + n, FILENAME_MAX, "xtradumps/cy%06lld", MetaData->CycleNumber);
//	n += snprintf(dirname + n, FILENAME_MAX, "xtradumps", MetaData->CycleNumber);
	snprintf(cmd, FILENAME_MAX + 9, "mkdir -p %s", dirname);
//	TRACEF("Creating directory '%s'...", dirname);
	system(cmd);

	snprintf(fname, FILENAME_MAX, "%s/%s-%06lld-%02lld-%lld-%02lld-%s.py", dirname, dumpPrefix, MetaData->CycleNumber,
				MyProcessorNumber, level, g->GetGridID(), dumpSuffix);
//	TRACEF("Creating file '%s'...", fname);
	file = fopen(fname, "w+");

	TRACEF("%p %p", file, g);

	fprintf(file, "%s\n"
			"#\n"
			"%s\n"
			"#\n"
			"# To load in python:  d = eval(open(path-to-this-file, 'r').read())\n"
			"# To examine the data, for example:  print(d.data[5:10], d.gridDims)\n"
			"#\n"
			"type('', (), dict(\n"
			"cycle = %lld,\n"
			"time = %f,\n"
			"cols = 'tag i j k'.split(),\n"
			"prefix = '%s',\n"
			"suffix = '%s',\n"
			"level = %lld,\n"
			"gridID = %lld,\n"
			"cellWidth = %e,\n"
			"nGhostZones = %lld,\n"
			"totalCells = %lld,\n"
			"gridDims = ( %lld , %lld , %lld ),\n"
			"gridLeftEdge = ( %e , %e , %e ),\n"
			"gridRightEdge = ( %e , %e , %e ),\n"
			"data = np.array([\n",
			dumpPreamble1, dumpPreamble2, MetaData->CycleNumber, MetaData->Time, dumpPrefix, dumpSuffix, level,
			g->GetGridID(), g->GetCellWidth(0, 0), NumberOfGhostZones, g->GetGridSize(), g->GetGridDimension(0),
			g->GetGridDimension(1), g->GetGridDimension(2), // in, jn, kn,
			g->GetGridDimension(0) * g->GetGridDimension(1) * g->GetGridDimension(2), // in, jn, kn,
			g->GetCellLeftEdge(0, 0), g->GetCellLeftEdge(1, 0), g->GetCellLeftEdge(2, 0),
			g->GetCellLeftEdge(0, g->GetGridDimension(0)), g->GetCellLeftEdge(1, g->GetGridDimension(1)),
			g->GetCellLeftEdge(2, g->GetGridDimension(2)));

	if(filename)
	{
		if(*filename == NULL)
			*filename = fname;
	}
	else
		delete fname;

	return file;
}

void negEFile_close(FILE *file, const char* filename)
{
	if(file && file != stderr && file != stdout)
	{
		fprintf(file, "])))\n");
//		TRACEF("Closing '%s'...", filename);
		fclose(file);
	}
}

int ClearOuterVelocities(float *u, float *v, float *w, float *totE, float *rhoField, float *pressure, int in, int jn,
int kn, int rank, float dx[], float dy[], float dz[], FLOAT** CellLeftEdge, int level, TopGridData *MetaData,
int gridID, grid *g, char* dumpPrefix, char *dumpSuffix, char* dumpPreamble)
{
//	arr_set(u, g->GetGridSize(), 0);
//	arr_set(v, g->GetGridSize(), 0);
//	arr_set(w, g->GetGridSize(), 0);
//	return SUCCESS;
#define DUMP_PREAMBLE2 "# The record format is" \
	"# tag, i, j, k\n" \
	"#   i,j,k are integer zone coordinates in the current grid;" \
	"#   tag = 1 -- total energy <= 0 before subtracting kinetic energy;\n" \
	"#   tag = 2 -- total energy <= 0 after subtracting kinetic energy;\n" \
	"#   tag = 0 -- N/A.\n"

	//	const double TINY_V = tiny_number;
	const FLOAT DX = dx[0];
	FLOAT r_x, r_y, r_z;
	size_t index;
	bool clearRadial, isOutsideR, clearTan, dampTan;
	double r, vdotr, rx, ry, rz, v_rx, v_ry, v_rz, v_tx, v_ty, v_tz;
	double oldVx, oldVy, oldVz, newVx, newVy, newVz;
	double oldKE, newKE;
	char negEFileName[FILENAME_MAX];
	char* negEptr = negEFileName;
	FILE *negEFile = NULL;
	bool negEFileCreate = (level >= 0) && MetaData && dumpPreamble && dumpPrefix && dumpSuffix;
	int negE1 = 0, negE2 = 0;
	int flagsIOT = 4 * OuterVelocitiesClearInward + 2 * OuterVelocitiesClearOutward + OuterVelocitiesClearTangential;
	if(OuterVelocitiesClearInward == 4 || OuterVelocitiesClearOutward == 4 || OuterVelocitiesClearTangential == 4)
		flagsIOT = 10;
	else if(OuterVelocitiesClearInward == 3 || OuterVelocitiesClearOutward == 3 || OuterVelocitiesClearTangential == 3)
		flagsIOT = 9;
	else if(OuterVelocitiesClearInward == 2 || OuterVelocitiesClearOutward == 2 || OuterVelocitiesClearTangential == 2)
		flagsIOT = 8;

	int edgeFilterZones = (OuterVelocitiesDistFromEdge > 0) ? nint(OuterVelocitiesDistFromEdge / DX) : 0;
	bool doAnyFilter, doSphereFilter, doEdgeFilter[6];

	TRACEF("  %p", MetaData);
	doSphereFilter = flagsIOT && (OuterVelocitiesSphereRadius > 0);
	doEdgeFilter[0] = (fabs(CellLeftEdge[0][NumberOfGhostZones] - DomainLeftEdge[0]) < DX / 8) && (edgeFilterZones > 0)
			&& (MetaData->LeftFaceBoundaryCondition[0] == 1);
	doEdgeFilter[2] = (fabs(CellLeftEdge[1][NumberOfGhostZones] - DomainLeftEdge[1]) < DX / 8) && (edgeFilterZones > 0)
			&& (MetaData->LeftFaceBoundaryCondition[1] == 1);
	doEdgeFilter[4] = (fabs(CellLeftEdge[2][NumberOfGhostZones] - DomainLeftEdge[2]) < DX / 8) && (edgeFilterZones > 0)
			&& (MetaData->LeftFaceBoundaryCondition[2] == 1);
	doEdgeFilter[1] = (fabs(CellLeftEdge[0][in - NumberOfGhostZones] - DomainRightEdge[0]) < DX / 8)
			&& (edgeFilterZones > 0) && (MetaData->RightFaceBoundaryCondition[0] == 1);
	doEdgeFilter[3] = (fabs(CellLeftEdge[1][jn - NumberOfGhostZones] - DomainRightEdge[1]) < DX / 8)
			&& (edgeFilterZones > 0) && (MetaData->RightFaceBoundaryCondition[1] == 1);
	doEdgeFilter[5] = (fabs(CellLeftEdge[2][kn - NumberOfGhostZones] - DomainRightEdge[2]) < DX / 8)
			&& (edgeFilterZones > 0) && (MetaData->RightFaceBoundaryCondition[2] == 1);

	doAnyFilter = doSphereFilter || doEdgeFilter[0] || doEdgeFilter[1] || doEdgeFilter[2] || doEdgeFilter[3]
			|| doEdgeFilter[4] || doEdgeFilter[5];

	if(!doAnyFilter)
		return SUCCESS;

	const size_t I1 = (doEdgeFilter[0]) ? (NumberOfGhostZones + edgeFilterZones) : 0;
	const size_t J1 = (doEdgeFilter[2]) ? (NumberOfGhostZones + edgeFilterZones) : 0;
	const size_t K1 = (doEdgeFilter[4]) ? (NumberOfGhostZones + edgeFilterZones) : 0;
	const size_t I2 = (doEdgeFilter[1]) ? (in - NumberOfGhostZones - edgeFilterZones) : in;
	const size_t J2 = (doEdgeFilter[3]) ? (jn - NumberOfGhostZones - edgeFilterZones) : jn;
	const size_t K2 = (doEdgeFilter[5]) ? (kn - NumberOfGhostZones - edgeFilterZones) : kn;

	v_rx = v_ry = v_rz = 0;
	//################################################################################
	// Filter velocity between the sphere and the edge filters (the sphere with radius
	// OuterVelocitiesSphereRadius and the box, which is the domain shrunk by edgeFilterZones.)
	for(int k = K1; k < K2; k++)
	{
		if(!doSphereFilter)
			break;

		for(int j = J1; j < J2; j++)
		{
			for(int i = I1; i < I2; i++)
			{
				r_x = CELLCENTER(0, i) - SphericalGravityCenter[0];
				r_y = CELLCENTER(1, j) - SphericalGravityCenter[1];
				r_z = CELLCENTER(2, k) - SphericalGravityCenter[2];

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
//					TRACEF("  %lld  %e  %e", isOutsideR , r , OuterVelocitiesSphereRadius);
				}

				// Do not change velocities inside the radius.
				if(!isOutsideR)
					continue;

				index = i + in * (j + jn * k);
				oldVx = u[index];
				oldVy = v[index];
				oldVz = w[index];
				oldKE = (square(oldVx) + square(oldVy) + square(oldVz)) / 2;
				if(flagsIOT == 8 || flagsIOT == 9 || flagsIOT == 10)
				{
					// Find the coordinates where the cell radius vector crosses the sphere.
					long ijk2[3];
					FLOAT crr = r / OuterVelocitiesSphereRadius;
					FLOAT r_xyz[3];
					r_xyz[0] = r_x / crr + SphericalGravityCenter[0];
					if(r_xyz[0] <= g->CellLeftEdge[0][0])
						continue;
					if(r_xyz[0] >= g->CellLeftEdge[0][g->GetGridDimension(0)])
						continue;
					r_xyz[1] = r_y / crr + SphericalGravityCenter[1];
					if(r_xyz[1] <= g->CellLeftEdge[1][0])
						continue;
					if(r_xyz[1] >= g->CellLeftEdge[1][g->GetGridDimension(1)])
						continue;
					r_xyz[2] = r_z / crr + SphericalGravityCenter[2];
					if(r_xyz[2] <= g->CellLeftEdge[2][0])
						continue;
					if(r_xyz[2] >= g->CellLeftEdge[2][g->GetGridDimension(2)])
						continue;

					// Find the cell on the sphere that is on the current radius vector.
					size_t index2 = g->get_ijk_index(ijk2, r_xyz);

					// use the cell values on the sphere index2] to calculate
					// new values of the current cell [index]:
					newVx = u[index2];
					newVy = v[index2];
					newVz = w[index2];
					if(flagsIOT == 9)
					{
						newVx /= crr;
						newVy /= crr;
						newVz /= crr;
					}
					newKE = (square(newVx) + square(newVy) + square(newVz)) / 2;
					u[index] = newVx;
					v[index] = newVy;
					w[index] = newVz;

					if(flagsIOT == 10)
					{
						if(rhoField)
							rhoField[index] = rhoField[index2];
						if(totE)
							totE[index] = totE[index2];
						if(pressure)
							pressure[index] = pressure[index2];
						newKE = oldKE = 0;
					}
				}
				else if(flagsIOT == 7)
				{
					// If clear the normal component regardless of direction
					// (i.e both inward and outward) and clear the tangential
					// then set velocity to zero and move to the next cell.
					newKE = u[index] = v[index] = w[index] = 0;
				}
				else
				{
					// Project the velocity onto the radius.
					vdotr = r_x * u[index] + r_y * v[index] + r_z * w[index];

					// Determine if the radial component should be cleared
					// based on the parameters and depending on whether it is
					// pointing toward the center (inward) or away (outward).
					clearRadial = (vdotr < 0) ? OuterVelocitiesClearInward : OuterVelocitiesClearOutward;

					clearTan = OuterVelocitiesClearTangential
							&& (OuterVelocitiesSphereRadius2 <= 0
									|| OuterVelocitiesSphereRadius2 <= OuterVelocitiesSphereRadius
									|| OuterVelocitiesSphereRadius2 <= r);
					dampTan = OuterVelocitiesClearTangential
							&& OuterVelocitiesSphereRadius2 > OuterVelocitiesSphereRadius
							&& OuterVelocitiesSphereRadius2 > 0 && OuterVelocitiesSphereRadius2 > r;

					if(clearTan && clearRadial)
					{
						// If clear both the tangential and the radial components,
						// set the velocity to zero and move to the next zone.
						newKE = u[index] = v[index] = w[index] = 0;
					}
					else if(!clearTan && !dampTan && !clearRadial)
					{
						// If clear neither tangential nor the radial component,
						// there is nothing to do for this zone.

						// Continue with the next zone, or ...
						continue;

						// ... or validate that the energy is positive at the end of the loop.
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
						FLOAT dampCoeff = 0;
						if(dampTan)
						{
							dampCoeff = square(
									M_PI_2 / (OuterVelocitiesSphereRadius2 - OuterVelocitiesSphereRadius)
											* (r - OuterVelocitiesSphereRadius));
							dampCoeff = 1 + dampCoeff * (-0.5 + square(dampCoeff) / 720); // approx Cosine.
						}
						r *= r;
						v_rx = vdotr * r_x / r;
						v_ry = vdotr * r_y / r;
						v_rz = vdotr * r_z / r;

						// In this branch we clear either the radial
						// or the tangential components of the velocity.
						// The neither-nor cases are handled in the above if/else-if.
						if(clearRadial)
						{
							newVx = oldVx - v_rx;
							newVy = oldVy - v_ry;
							newVz = oldVz - v_rz;
							if(dampTan)
							{
								newVx *= dampCoeff;
								newVy *= dampCoeff;
								newVz *= dampCoeff;
							}
						}
						else if(clearTan)
						{
							newVx = v_rx;
							newVy = v_ry;
							newVz = v_rz;
						}
						else
						{
							newVx = (oldVx - v_rx) * dampCoeff + v_rx;
							newVy = (oldVy - v_ry) * dampCoeff + v_ry;
							newVz = (oldVz - v_rz) * dampCoeff + v_rz;
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
				int fixType = 0; // Nothing to fix
				if(oldE <= 0)
				{
					newE = (newKE > 0) ? newKE : (oldKE > 0) ? oldKE : 0 + tiny_pressure / (1 - Gamma);
					fixType = 1;
					negE1++;
				}
				else if(newE <= 0)
				{
					newE = tiny_pressure / (1 - Gamma);
					if(newE < oldE)
						newE = oldE;
					fixType = 2;
					negE2++;
				}
				totE[index] = newE;
				if(fixType && negEFileCreate)
				{
					if(negEFile == NULL)
						negEFile = negEFile_open(&negEptr, MetaData, level, g, "ne", dumpSuffix, dumpPreamble,
						DUMP_PREAMBLE2);
					fprintf(negEFile, "(%lld,%lld,%lld,%lld),\n", fixType, i, j, k);
				}
			}
		}
	}

	//################################################################################
	// Process zones for each of the edge filters present.
	for(int edge = 0; edge < 6; edge++)
	{
		if(!doEdgeFilter[edge])
			continue;

		const bool isRight = edge % 2;
		const int edgeDim = edge / 2;
		const size_t I1b = (edge == 1) ? I2 : 0;
		const size_t J1b = (edge == 3) ? J2 : 0;
		const size_t K1b = (edge == 5) ? K2 : 0;
		const size_t I2b = (edge == 0) ? I1 : in;
		const size_t J2b = (edge == 2) ? J1 : jn;
		const size_t K2b = (edge == 4) ? K1 : kn;
		const int sourceLayer = (edge == 0) ? I2b : (edge == 1) ? (I1b - 1) : (edge == 2) ? J2b :
								(edge == 3) ? (J1b - 1) : (edge == 4) ? K2b : (K1b - 1);

		v_rx = v_ry = v_rz = 0;
		// Clear velocity beyond a certain radius
		for(int k = K1b; k < K2b; k++)
		{
			for(int j = J1b; j < J2b; j++)
			{
				for(int i = I1b; i < I2b; i++)
				{
					FLOAT r_xyz[3], r_xyz2[3];
					r_xyz[0] = CELLCENTER(0, i) - SphericalGravityCenter[0];
					r_xyz[1] = CELLCENTER(1, j) - SphericalGravityCenter[1];
					r_xyz[2] = CELLCENTER(2, k) - SphericalGravityCenter[2];
					arr_set(r_xyz2, 3, 0);
					FLOAT c = r_xyz2[edgeDim] = CELLCENTER(edgeDim, sourceLayer) - SphericalGravityCenter[edgeDim];
					c /= r_xyz[edgeDim];
					for(int dim = 0; dim < 3; dim++)
						if(dim != edgeDim)
							r_xyz2[dim] = r_xyz[dim] * c;

					index = i + in * (j + jn * k);
					long ijk2[3];
					size_t index2 = g->get_ijk_index(ijk2, r_xyz2);

					oldVx = u[index];
					oldVy = v[index];
					oldVz = w[index];
					oldKE = lensquaredl(oldVx, oldVy, oldVz) / 2;

					newVx = u[index2];
					newVy = v[index2];
					newVz = w[index2];
					newVx *= c * c;
					newVy *= c * c;
					newVz *= c * c;
					newKE = lensquaredl(newVx, newVy, newVz) / 2;
					u[index] = newVx;
					v[index] = newVy;
					w[index] = newVz;

					float dE = newKE - oldKE;
					float oldE = totE[index];
					float newE = oldE + dE;
					int fixType = 0; // Nothing to fix
					if(oldE <= 0)
					{
						newE = (newKE > 0) ? newKE : (oldKE > 0) ? oldKE : 0 + tiny_pressure / (1 - Gamma);
						fixType = 1;
						negE1++;
					}
					else if(newE <= 0)
					{
						newE = tiny_pressure / (1 - Gamma);
						if(newE < oldE)
							newE = oldE;
						fixType = 2;
						negE2++;
					}
					totE[index] = newE;
					if(fixType && negEFileCreate)
					{
						if(negEFile == NULL)
							negEFile = negEFile_open(&negEptr, MetaData, level, g, "ne", dumpSuffix, dumpPreamble,
							DUMP_PREAMBLE2);
						fprintf(negEFile, "(%lld,%lld,%lld,%lld),\n", fixType, i, j, k);
					}
				}
			}
		}
	}

	if(negE1 + negE2)
	{
		TRACEF(" %s: Bad energy in %lld zones of total %lld (%lld before, %lld after correction)", negEFileName,
				negE1 + negE1 + negE2, g->GetGridSize(), negE1, negE2);
	}

	if(negEFileCreate)
		negEFile_close(negEFile, negEFileName);

	return SUCCESS;
#undef DUMP_PREAMBLE2
}

int grid::ClearOuterVelocities(float *pressure, int level, TopGridData *MetaData, char* dumpPrefix, char *dumpSuffix,
	char* dumpPreamble)
{
	if(MyProcessorNumber != ProcessorNumber)
		return SUCCESS;

	float *rhoField, *vxField, *vyField, *vzField, *totEField;
	MHD_SNIA_GetFields(&rhoField, &totEField, NULL, &vxField, &vyField, &vzField, NULL, NULL, NULL, NULL, NULL,
	NULL,
						NULL,
						NULL);

	return ::ClearOuterVelocities(vxField, vyField, vzField, totEField, rhoField, pressure, GridDimension[0],
									GridDimension[1], GridDimension[2], GridRank, CellWidth[0], CellWidth[1],
									CellWidth[2], CellLeftEdge, level, MetaData, this->ID, this, dumpPrefix, dumpSuffix,
									dumpPreamble);
}

int ClearOuterVelocities(float *u, float *v, float *w, float *totE, int in, int jn, int kn, int rank,
float dx[],
float dy[], float dz[], FLOAT** CellLeftEdges, int level, TopGridData *MetaData, int gridID, grid *g, char* dumpPrefix,
	char *dumpSuffix, char* dumpPreamble)
{
	return ClearOuterVelocities(u, v, w, totE, NULL, NULL, in, jn, kn, rank, dx, dy, dz, CellLeftEdges, level, MetaData,
								gridID, g, dumpPrefix, dumpSuffix, dumpPreamble);
}
