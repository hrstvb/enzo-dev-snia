/***********************************************************************
 /
 /  GRID CLASS (WRAPPER FOR ZEUS HYDRO)
 /
 /  written by: Greg Bryan
 /  date:       November, 1994
 /  modified1:
 /
 /  PURPOSE:
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 c
 c  INPUTS:
 c     d       - density field (includes boundary zones)
 c     dx,y,z  - zone width arrays for each dimension
 c     e       - total specific energy field
 c     end     - array (of dimension 3) specifying the end of the active
 c               region for reach dimension (zero based)
 c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
 c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
 c     ge      - gas energy (used when idual = 1)
 c     gr_x,y,zacc - gravitational acceleration fields
 c     gravity - flag indicating whether or not to use gravity field (1 = yes)
 c     gridvel - bulk grid velocity (vector of length 3)
 c     i,j,kn  - dimensions of field arrays
 c     idiff   - diffusion flag (0 = off)
 c     idual   - dual energy formalism flag (0 = off)
 c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
 c     nhy     - cycle number (for better operator splitting)
 c     rank    - dimension of problem (not currently used)
 c     start   - array (of dimension 3) specifying the start of the active
 c               region fo reach dimension (zero based)
 c     tmp     - temporary work space (30 * largest_slice)
 c     u       - x-velocity field
 c     v       - y-velocity field
 c     w       - z-velocity field
 c     bottom  - true (1) if this is the lowest level
 c     minsupecoef - coefficient for minimum pressure support (0 - not used)
 /
 ************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//
#include <stdio.h>
#include <math.h>
#include "DebugMacros.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

int Zeus_xTransport(float *d, float *e, float *u, float *v, float *w,
int in, int jn, int kn, int rank,
int is, int ie, int js, int je, int ks, int ke,
float dt, float dx[], float *f1, int bottom,
int nsubgrids, long_int GridGlobalStart[], fluxes *SubgridFluxes[], int DensNum, int TENum,
int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
int NumberOfColours, int colnum[]);

int Zeus_yTransport(float *d, float *e, float *u, float *v, float *w,
int in, int jn, int kn, int rank,
int is, int ie, int js, int je, int ks, int ke,
float dt, float dy[], float *f1, int bottom,
int nsubgrids, long_int GridGlobalStart[], fluxes *SubgridFluxes[], int DensNum, int TENum,
int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
int NumberOfColours, int colnum[]);

int Zeus_zTransport(float *d, float *e, float *u, float *v, float *w,
int in, int jn, int kn, int rank,
int is, int ie, int js, int je, int ks, int ke,
float dt, float dz[], float *f1, int bottom,
int nsubgrids, long_int GridGlobalStart[], fluxes *SubgridFluxes[], int DensNum, int TENum,
int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
int NumberOfColours, int colnum[]);

bool debugprintboundaryForLevel(int level)
{
	return level != 0;
}

int debugprintboundary(char* LABEL, char* label, float *f, FLOAT* F, char* label2, float *f2, FLOAT* F2, int in, int jn,
int kn,
int rank, int is, int ie, int js,
int je, int ks, int ke, float dt, float dx[], float dy[], float dz[], FLOAT** CellLeftEdges, int level, int gridId,
float time, int cycle)
{
	if(debugprintboundaryForLevel(level))
		return 0;

	int IB = 0;
	int IE = 4;
	int JB = 0;
	int JE = 5;
	int KB = 0;
	int KE = 5;
#define MARKER "MARKER|"
	fprintf(stderr, MARKER "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fprintf(stderr, MARKER "+++++ level=%lld, gridID=%lld, t=%7.5f  :+++++++++++++\n", level, gridId, time);
	for(size_t i = IB; i <= IE; i++)
	{
		if(label2)
			fprintf(stderr,
					MARKER "  i=%lld :----------  %s  :----  level=%lld  %s  :------------                        %s  :--------------------------\n",
					i, label, level, LABEL, label2);
		else
			fprintf(stderr, MARKER "  i=%lld :----------  %s  :----  level=%lld  %s\n", i, label, level, LABEL);
		for(size_t j = JB; j <= JE; j++)
		{
			fprintf(stderr, MARKER "  j=%lld, k=%lld..%lld:", j, KB, KE);
			for(size_t k = KB; k <= KE; k++)
			{
				size_t index = i + in * (j + jn * k);
				if(f)
					fprintf(stderr, "%11.3e", f[index]);
				else if(F)
					fprintf(stderr, "  %e", F[index]);
			}
			if(f2 || F2)
			{
				fprintf(stderr, "            ");
				for(size_t k = KB; k <= KE; k++)
				{
					size_t index = i + in * (j + jn * k);
					if(f2)
						fprintf(stderr, "%11.3e", f2[index]);
					else if(F2)
						fprintf(stderr, "  %e", F2[index]);
				}
			}
			fprintf(stderr, "\n");
		}
	}
}

int debugprintboundary(char* LABEL, char* label, float *f, FLOAT* F, int in, int jn, int kn,
int rank, int is, int ie, int js,
int je, int ks, int ke, float dt, float dx[], float dy[], float dz[], FLOAT** CellLeftEdges, int level, int gridId,
float time, int cycle)
{
	return debugprintboundary(LABEL, label, f, F, NULL, NULL, NULL, in, jn, kn, rank, is, ie, js, je, ks, ke, dt, dx,
								dy, dz, CellLeftEdges, level, gridId, time, cycle);
}

int debugprintboundary(char* LABEL, float *d, float *e, float *u, float *v, float *w, float *p, int in, int jn, int kn,
int rank, int is, int ie, int js, int je, int ks, int ke, float dt, float dx[], float dy[], float dz[],
FLOAT** CellLeftEdges, int level, int gridId, float time, int cycle, float *gr_xacc, float *gr_yacc, float *gr_zacc)
{
	if(debugprintboundaryForLevel(level))
		return 0;

	fprintf(stderr, "BEGIN BOUNDARY %s === level=%lld, gridID=%lld, t=%7.5f ===========\n", LABEL, level, gridId, time);
	debugprintboundary(LABEL, "density", d, NULL, "u", u, NULL, in, jn, kn, rank, is, ie, js, je, ks, ke, dt, dx, dy,
						dz, CellLeftEdges, level, gridId, time, 0);
	debugprintboundary(LABEL, "v", v, NULL, "w", w, NULL, in, jn, kn, rank, is, ie, js, je, ks, ke, dt, dx, dy, dz,
						CellLeftEdges, level, gridId, time, 0);
	debugprintboundary(LABEL, "energy", e, NULL, in, jn, kn, rank, is, ie, js, je, ks, ke, dt, dx, dy, dz,
						CellLeftEdges, level, gridId, time, 0);
	debugprintboundary(LABEL, "pressure", p, NULL, "gravity_x", gr_xacc, NULL, in, jn, kn, rank, is, ie, js, je, ks, ke,
						dt, dx, dy, dz, CellLeftEdges, level, gridId, time, 0);
	debugprintboundary(LABEL, "gravity_y", gr_yacc, NULL, "gravity_z", gr_zacc, NULL, in, jn, kn, rank, is, ie, js, je,
						ks, ke, dt, dx, dy, dz, CellLeftEdges, level, gridId, time, 0);
	fprintf(stderr, "END BOUNDARY %s ===========================================\n", LABEL);

	return 0;
}

int ZeusSource(float *d, float *e, float *u, float *v, float *w, float *p, float *cr,
int in, int jn, int kn, int rank, int igamfield,
int is, int ie, int js, int je, int ks, int ke,
float C1, float C2, int ipresfree,
float *gamma, float dt, float pmin, float dx[], float dy[], float dz[], FLOAT** CellLeftEdges,
int gravity, float *gr_xacc, float *gr_yacc, float *gr_zacc,
int bottom, float minsupecoef, int CRModel, float CRgamma, size_t index8, size_t i8, size_t j8, size_t k8,
	size_t index9, size_t i9, size_t j9, size_t k9,
	FLOAT x8,
	FLOAT y8,
	FLOAT z8,
	FLOAT x9,
	FLOAT y9,
	FLOAT z9,
	int level, int gridId, float time, int cycle, int computePressureOnly);

int GetUnits(float *DensityUnits, float *LengthUnits,
float *TemperatureUnits, float *TimeUnits,
float *VelocityUnits, double *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::ZeusSolver(float *gamma, int igamfield, int nhy,
float dx[], float dy[], float dz[],
int gravity, int NumberOfSubgrids, long_int GridGlobalStart[], fluxes *SubgridFluxes[],
int NumberOfColours, int colnum[], int bottom,
float minsupecoef)
{

	/*  Locals */

	int i, ie, is, j, je, js, k, ks, ke, n, ixyz, ret;
	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum;
	float pmin;
	float *d, *e, *u, *v, *w, *cr, *m;

	/* Error Check */

	for(i = 0; i < GridRank; i++)
		if(GridDimension[i] > MAX_ANY_SINGLE_DIRECTION)
		{
			ENZO_FAIL("ZEUS_MAIN: A grid dimension is too long (increase max_any_single_direction.)\n");
		}

	/* Allocate temporary space for Zeus_Hydro. */

	int size = GridDimension[0] * GridDimension[1] * GridDimension[2];
	float *p = new float[size];

	/* Find fields: density, total energy, velocity1-3 and set pointers to them
	 Create zero fields for velocity2-3 for low-dimension runs because solver
	 assumes they exist. */

	this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);
	if(CRModel)
	{
		if((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
			ENZO_FAIL("Cannot Find Cosmic Rays");
		cr = BaryonField[CRNum];
	}

	d = BaryonField[DensNum];
	e = BaryonField[TENum];
	u = BaryonField[Vel1Num];
	v = BaryonField[Vel2Num];
	w = BaryonField[Vel3Num];
	if(GridRank < 2)
	{
		v = new float[size];
		for(i = 0; i < size; i++)
			v[i] = 0;
	}
	if(GridRank < 3)
	{
		w = new float[size];
		for(i = 0; i < size; i++)
			w[i] = 0;
	}

	/*  Set grid start and end indicies */

	is = GridStartIndex[0];
	js = GridStartIndex[1];
	ks = GridStartIndex[2];
	ie = GridEndIndex[0];
	je = GridEndIndex[1];
	ke = GridEndIndex[2];

	/***********************************************/
	/***********************************************/
	FLOAT DX = dx[0];
	int level = nint(log2(TopGridDx[0] / DX));
	FLOAT xyz8[3] = { -600e5 + DX / 2, -600e5 + DX / 2, -600e5 + DX / 2 };
	FLOAT &x8 = xyz8[0];
	FLOAT &y8 = xyz8[1];
	FLOAT &z8 = xyz8[2];
	long ijk8[3];
	long &i8 = ijk8[0];
	long &j8 = ijk8[1];
	long &k8 = ijk8[2];
	size_t index8 = get_ijk_index(ijk8, xyz8);
	if(!PointInGridActiveNB(xyz8))
		index8 = 0;

	FLOAT xyz9[3] = { 600e5 - DX / 2, 600e5 - DX / 2, 600e5 - DX / 2 };
	FLOAT &x9 = xyz9[0];
	FLOAT &y9 = xyz9[1];
	FLOAT &z9 = xyz9[2];
	long ijk9[3];
	long & i9 = ijk9[0];
	long & j9 = ijk9[1];
	long & k9 = ijk9[2];
	size_t index9 = get_ijk_index(ijk9, xyz9);
	if(!PointInGridActiveNB(xyz9))
		index9 = 0;

	/***********************************************/
	/***********************************************/

	//  If NumberOfGhostZones is set to 4, then use the extra space
	if(is == 4)
	{
		is = is - 1;
		ie = ie + 1;
	}
	if(js == 4)
	{
		js = js - 1;
		je = je + 1;
	}
	if(ks == 4)
	{
		ks = ks - 1;
		ke = ke + 1;
	}

	/* Set minimum pressure (better if it were a parameter) */

	pmin = tiny_number;

	/* Error check */

	for(i = 0; i < size; i++)
	{
		if(fabs(u[i]) > dx[0] / dtFixed || fabs(v[i]) > dy[0] / dtFixed || fabs(w[i]) > dz[0] / dtFixed)
		{
			fprintf(stderr, "u,v,w,d,e=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM"  dt=%"GSYM"\n", u[i], v[i],
					w[i], d[i], e[i], dx[0], dtFixed);
			ENZO_FAIL("Velocity too fast! (pre-call)\n");
		}
	}

	ZeusSource(d, e, u, v, w, p, cr, GridDimension[0], GridDimension[1], GridDimension[2], GridRank, igamfield, is, ie,
				js, je, ks, ke, ZEUSLinearArtificialViscosity, ZEUSQuadraticArtificialViscosity, PressureFree, gamma,
				dtFixed, pmin, dx, dy, dz, CellLeftEdge, gravity, AccelerationField[0], AccelerationField[1],
				AccelerationField[2], bottom, minsupecoef, CRModel, CRgamma, index8, i8, j8, k8, index9, i9, j9, k9,
				x8, y8, z8, x9, y9, z9, level, ID, dtFixed, 0, 1);

	debugprintboundary("BEFORE SOURCE", d, e, u, v, w, p, GridDimension[0], GridDimension[1], GridDimension[2],
						GridRank, is, ie, js, je, ks, ke, dtFixed, dx, dy, dz, CellLeftEdge, level, ID, dtFixed, 0,
						AccelerationField[0], AccelerationField[1], AccelerationField[2]);

	/*   1) Add source terms */

	/***********************************************/
	/***********************************************/
	if(ZeusSource(d, e, u, v, w, p, cr, GridDimension[0], GridDimension[1], GridDimension[2], GridRank, igamfield, is,
					ie, js, je, ks, ke, ZEUSLinearArtificialViscosity, ZEUSQuadraticArtificialViscosity, PressureFree,
					gamma, dtFixed, pmin, dx, dy, dz, CellLeftEdge, gravity, AccelerationField[0], AccelerationField[1],
					AccelerationField[2], bottom, minsupecoef, CRModel, CRgamma, index8, i8, j8, k8, index9, i9, j9,
					k9, x8, y8, z8, x9, y9, z9, level, ID, dtFixed, 0, 0) == FAIL)
	{
		fprintf(stderr, "P(%"ISYM"): Error in ZeusSource on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber, nhy,
				dtFixed);
		fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1],
				GridDimension[2]);
		ENZO_FAIL("Error in ZeusSource!\n");
	}

	debugprintboundary("AFTER SOURCE", d, e, u, v, w, p, GridDimension[0], GridDimension[1], GridDimension[2], GridRank,
						is, ie, js, je, ks, ke, dtFixed, dx, dy, dz, CellLeftEdge, level, ID, dtFixed, 0,
						AccelerationField[0], AccelerationField[1], AccelerationField[2]);

	/* Error check */

	float CRcs = 0.0;
	if(CRmaxSoundSpeed != 0.0)
	{
		// Get system of units
		float CRsound, DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, Time;
		double MassUnits;
		if(GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits,
					Time) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		CRsound = CRmaxSoundSpeed / VelocityUnits;
		CRcs = (CRgamma - 1.0) / (CRsound * CRsound);
	}

	for(i = 0; i < size; i++)
	{
		if(fabs(u[i]) > dx[0] / dtFixed || fabs(v[i]) > dy[0] / dtFixed || fabs(w[i]) > dz[0] / dtFixed)
		{
			fprintf(stderr, "u,v,w,d,e=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM"  dt=%"GSYM"\n", u[i], v[i],
					w[i], d[i], e[i], dx[0], dtFixed);
			ENZO_FAIL("Velocity too fast! (post-call)\n");
		}

		/* -- density/TE floor for CR model -- */

		if(CRModel)
		{
			if(CRdensFloor != 0.0 && d[i] < CRdensFloor)
				d[i] = CRdensFloor;
			if(CRcs != 0.0 && d[i] < CRcs * cr[i])
				d[i] = CRcs * cr[i];   // Limits sound-speed
			if(e[i] < tiny_number * 1e-5)
				e[i] = tiny_number * 1e-5;
		} // end cr model if
	} // end i for

	/*  2) Transport step */

	ixyz = nhy % GridRank;
	for(n = ixyz; n <= ixyz + GridRank - 1; n++)
	{

		/* Transport step - x direction */

		if((n % GridRank) == 0)
			ret = Zeus_xTransport(d, e, u, v, w, GridDimension[0], GridDimension[1], GridDimension[2], GridRank, is, ie,
									js, je, ks, ke, dtFixed, dx, p, bottom, NumberOfSubgrids, GridGlobalStart,
									SubgridFluxes, DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, BaryonField,
									NumberOfColours, colnum);

		/*  Transport step - y direction */

		if((n % GridRank) == 1 && GridRank > 1)
			ret = Zeus_yTransport(d, e, u, v, w, GridDimension[0], GridDimension[1], GridDimension[2], GridRank, is, ie,
									js, je, ks, ke, dtFixed, dy, p, bottom, NumberOfSubgrids, GridGlobalStart,
									SubgridFluxes, DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, BaryonField,
									NumberOfColours, colnum);

		/*  Transport step - z direction */

		if((n % GridRank) == 2 && GridRank > 2)
			ret = Zeus_zTransport(d, e, u, v, w, GridDimension[0], GridDimension[1], GridDimension[2], GridRank, is, ie,
									js, je, ks, ke, dtFixed, dz, p, bottom, NumberOfSubgrids, GridGlobalStart,
									SubgridFluxes, DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, BaryonField,
									NumberOfColours, colnum);

		if(ret == FAIL)
		{
			fprintf(stderr, "P(%"ISYM"): Error on ZeusTransport dim=%"ISYM" (Cycle = %"ISYM", dt=%"GSYM")\n",
					MyProcessorNumber, n % GridRank, nhy, dtFixed);
			fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1],
					GridDimension[2]);
			ENZO_FAIL("Error in ZeusSource!\n");
		}

	} // end loop over n

	/* Clean up */
	delete[] p;
	if(GridRank < 2)
		delete[] v;
	if(GridRank < 3)
		delete[] w;

	return SUCCESS;

}
