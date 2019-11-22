/***********************************************************************
 /
 /  GRID CLASS (INITIALIZE THE MHD BLAST GRID)
 /
 /  written by: David Collins
 /  date:       2004-2013
 /
 /  modified1:  Boyan Hristov
 /  date:       2015
 /  		Added perturbation method 33 creating a plane interface
 /		between zones A and B perturbed by a sine offset.
 /
 /  PURPOSE:  See MHDBlastInitialize.C for parameters
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <cctype> //isspace
#include <stdlib.h>
#include <cmath>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "mylimiters.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

#include "DebugMacros.h"
#include "MHDInitialProfile.h"
#include "TriSpherePerturb.C"

//#ifndef SRC_ENZO_INITIALPROFILE_H_
//#define SRC_ENZO_INITIALPROFILE_H_
//
//#define LINE_MAX_LENGTH (512)
//#define CMP_A_B(a, b) ((b>a)-(b<a))
//#define SIGN_A(a) (CMP_A_B(0,a))
//
//#define PROFILE_ASC_SORT (1)
//#define PROFILE_DESC_SORT (-1)
//#define PROFILE_UNSORTED (0)
//
//#define PROFILE_EMPTY_LINE (0)
//#define PROFILE_COMMENT_LINE (1)
//#define PROFILE_COL_NAMES_LINE (2)
//#define PROFILE_DATA_LINE (3)
//#define PROFILE_TIME_LINE (4)

//using namespace std;

int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z);
float SphericalGravityGetAt(FLOAT r);
//void WriteInitialProfile(char* name, FLOAT* RR, FLOAT* RHO, FLOAT* GG, FLOAT* PP, FLOAT* UU, size_t n, FLOAT K,
//FLOAT gamma);

inline float internalEnergy(float rho, float rhoNi, MHDInitialProfile* profile, FLOAT radius)
{
	double T = 0, E = 0;
	if(profile->interpolateTemperature(&T, radius))
		T = (rhoNi / rho > 0.5) ? 5e9 : 1e6;
	E += 1.5 * (1 - 0.75 * rhoNi / rho) * R_gas * T / 14;
	E += EOSPolytropicFactor * pow(rho, Gamma - 1) / (Gamma - 1);
	return E;
}

int polytropicPressureAtSmallR(float* P, float* rho, FLOAT r, float rho_c, float * P_c = NULL, double* ksi = NULL,
	double* ksiFactor = NULL)
{
	float P_c2, P2, rho2;
	double ksi2, ksi4, ksiFactor2;
	int err;

	// Initialize the center cells
	double KPolytropic = EOSPolytropicFactor;
	double nPolytropic = 1 / (Gamma - 1);

	// Find the dimensionless radius for the expansions further below:
	ksiFactor2 = (nPolytropic + 1) * KPolytropic * POW(rho_c, 1 / nPolytropic - 1);
	ksiFactor2 /= 4 * M_PI * SphericalGravityConstant;
	ksiFactor2 = sqrt(ksiFactor2);
	ksi2 = r / ksiFactor2;

	if(ksi)
		*ksi = ksi2;
	if(ksiFactor)
		*ksiFactor = ksiFactor2;

	ksi2 *= ksi2;
	ksi4 = ksi2 * ksi2;

	if(rho)
	{
		rho2 = 1;
		rho2 += -nPolytropic / 6 * ksi2;
		rho2 += -nPolytropic * (22 * nPolytropic + 5) / 360 * ksi4;
		rho2 *= rho_c;

		*rho = rho2;
	}

	if(P || P_c)
	{
		P_c2 = KPolytropic * POW(rho_c, Gamma);
		if(P_c)
			*P_c = P_c2;

		if(P)
		{	// Expand the pressure for small radius:
//			P2 = 1; // 0th order
//			P2 += -(nPolytropic + 1) / 6 * ksi2; // 2nd order
//			P2 += -11.0 * nPolytropic * (nPolytropic + 1) / 180.0 * ksi4; // 4th order
			double P2_0 = 1; // 0th order
			double P2_2 = -(nPolytropic + 1) / 6 * ksi2; // 2nd order
			double P2_4 = -11.0 * nPolytropic * (nPolytropic + 1) / 180.0 * ksi4; // 4th order
			P2 = P2_0 + P2_2 + P2_4;
//			TRACEF("Polytropic pressure ksi^2=%e, P0+P2+P4 = %e + %e + %e = %e, P_c = %e, P2*P_c=%e, rho_c=%e", ksi2,
//					P2_0, P2_2, P2_4, P2, P_c2, P2 * P_c2, rho_c);
			P2 *= P_c2;

			*P = P2;
		}
	}

//	TRACEF("ijk=%lld %lld %lld xyz=%e %e %e r,xi,r/xi=%e %e %e   P=%e   P/P_c=%e", i, j, k, x, y, z, r, ksi, ksifactor,
//			P0, P1);

	return 0;
}

int polytropicPressureAtSmallR(float* P, float* rho, FLOAT x, FLOAT y, FLOAT z, double rho_c, double* P_c = NULL,
	double* ksi = NULL, double* ksiFactor = NULL, FLOAT *r = NULL)
{
	x -= SphericalGravityCenter[0];
	y -= SphericalGravityCenter[1];
	z -= SphericalGravityCenter[2];
	FLOAT r2 = lenl(x, y, z);
	if(r)
		*r = r2;

	return polytropicPressureAtSmallR(P, rho, r2, rho_c, P_c, ksi, ksiFactor);
}

int polytropicPressureAtSmallR(float* P, float* rho, int i, int j, int k, grid* g, double rho_c, double* P_c = NULL,
FLOAT* ksi = NULL, FLOAT* ksiFactor = NULL, FLOAT *r = NULL)
{
	FLOAT x = g->GetCellCenter(0, i);
	FLOAT y = g->GetCellCenter(1, j);
	FLOAT z = g->GetCellCenter(2, k);

	return polytropicPressureAtSmallR(P, rho, x, y, z, rho_c, P_c, ksi, ksiFactor, r);
}

double zeusPressure(double P_prev, double g, double dx, double rho_prev, double rho_this)
{
	double dP, P_result;
	dP = g * dx * (rho_prev + rho_this) / 2;
	P_result = P_prev + dP;
	return P_result;
}

double zeusPressure(size_t i_dest, size_t j_dest, size_t k_dest, int dir, int zeusDim, double dx, float* P_FIELD,
float* RHO_FIELD, grid* GRID, bool debug, int level, char* label)
{
	int gridId = GRID->GetGridID();
	size_t i_prev = i_dest - dir * (zeusDim == 0);
	size_t j_prev = j_dest - dir * (zeusDim == 1);
	size_t k_prev = k_dest - dir * (zeusDim == 2);
	size_t index_dest = GRID->GetIndex(i_dest, j_dest, k_dest);
	size_t index_prev = GRID->GetIndex(i_prev, j_prev, k_prev);

	double P_result;
	double P_prev = P_FIELD[index_prev];
	double rho_dest = RHO_FIELD[index_dest];
	double rho_prev = RHO_FIELD[index_prev];

	double x = (zeusDim == 0) ? GRID->GetCellLeftEdge(0, (dir > 0) ? i_dest : i_prev) : GRID->GetCellCenter(0, i_dest);
	double y = (zeusDim == 1) ? GRID->GetCellLeftEdge(1, (dir > 0) ? j_dest : j_prev) : GRID->GetCellCenter(1, j_dest);
	double z = (zeusDim == 2) ? GRID->GetCellLeftEdge(2, (dir > 0) ? k_dest : k_prev) : GRID->GetCellCenter(2, k_dest);
	x -= SphericalGravityCenter[0];
	y -= SphericalGravityCenter[1];
	z -= SphericalGravityCenter[2];
	double r = lenl(x, y, z);
	double r_dim = (zeusDim == 0) ? x : ((zeusDim == 1) ? y : z);
	double g = SphericalGravityGetAt(r);
	double g_dim = -r_dim / r * g;

	P_result = zeusPressure(P_prev, g_dim, dir * dx, rho_prev, rho_dest);
//						if(j == 8 && k == 8)
//						if(err)
//						{
//							TRACEF("ijk=%lld %lld %lld, xyzr=%e %e %e %e, g,dx=%e %e, P1,2=%e %e  d1,2=%e %e", i, j, k,
//									x, y, z, r, g, dx, P1, P2, rho1, rho2);
//						}

	debug &= i_dest <= 4 && (j_dest == 4) && (k_dest == 4);
	if(debug)
	{
		TRACEF("INITINITINIT  %s  %lld:%lld  %lld    %lld %lld %lld    %lld %lld %lld    P=%e  P1=%e  dP=%e    d=%e  d1=%e  da=%e    h=%e    dP/da/h=%e",
				label, level, gridId, zeusDim, i_dest, j_dest, k_dest, i_prev, j_prev, k_prev, P_result, P_prev,
				P_result - P_prev, rho_dest, rho_prev, 0.5 * (rho_dest + rho_prev), dx,
				(P_result - P_prev) / (dx * 0.5 * (rho_dest + rho_prev)));
		TRACEF("INITINITINIT  %s  %lld:%lld  %lld    %lld %lld %lld    %lld %lld %lld    g( %e, %e, %e )=%e    gimd=-g*( %e / %e )=%e",
				label, level, gridId, zeusDim, i_dest, j_dest, k_dest, i_prev, j_prev, k_prev, x, y, z, g, r_dim, r,
				g_dim);
	}

	return P_result;
}

double rkPressure(double P_prev2, double g_prev, double dx, double rho_prev)
{
	double dP, P_result;
	dP = 2 * g_prev * dx * rho_prev;
	P_result = P_prev2 + dP;
	return P_result;
}

double rkPressure(size_t i_dest, size_t j_dest, size_t k_dest, int integrationDirection, int integrationDimension,
	double dx, float* P_FIELD, float* RHO_FIELD, grid* GRID, bool debug, int level, char* label)
{
	int gridId = GRID->GetGridID();
	size_t i_prev2 = i_dest - 2 * integrationDirection * (integrationDimension == 0);
	size_t j_prev2 = j_dest - 2 * integrationDirection * (integrationDimension == 1);
	size_t k_prev2 = k_dest - 2 * integrationDirection * (integrationDimension == 2);
	size_t i_prev = i_dest - integrationDirection * (integrationDimension == 0);
	size_t j_prev = j_dest - integrationDirection * (integrationDimension == 1);
	size_t k_prev = k_dest - integrationDirection * (integrationDimension == 2);
	size_t index_dest = GRID->GetIndex(i_dest, j_dest, k_dest);
	size_t index_prev = GRID->GetIndex(i_prev, j_prev, k_prev);
	size_t index_prev2 = GRID->GetIndex(i_prev2, j_prev2, k_prev2);

	double P_result;
	double rho_prev2 = RHO_FIELD[index_prev2];
	double P_prev2 = rho_prev2 * P_FIELD[index_prev2];
	double rho_prev = RHO_FIELD[index_prev];

	double x = GRID->GetCellCenter(0, i_prev);
	double y = GRID->GetCellCenter(1, j_prev);
	double z = GRID->GetCellCenter(2, k_prev);
	x -= SphericalGravityCenter[0];
	y -= SphericalGravityCenter[1];
	z -= SphericalGravityCenter[2];
	double r = lenl(x, y, z);
	double r_dim = (integrationDimension == 0) ? x : ((integrationDimension == 1) ? y : z);
	double g = SphericalGravityGetAt(r);
	double g_dim = -r_dim / r * g;

	P_result = rkPressure(P_prev2, g_dim, integrationDirection * dx, rho_prev);

	if(debug)
	{
		TRACEF("ijk=%lld %lld %lld   xyz=%e %e %e   " "r,r_dim,=%e %e   g,g_dim=%e %e   P,P_prev2,rho=%e %e %e" "   ijk_prev=%lld %lld %lld" "   ijk_prev2=%lld %lld %lld    %lld %e",
				i_dest, j_dest, k_dest, x, y, z, r, r_dim, g, g_dim, P_result, P_prev2, rho_prev, i_prev, j_prev,
				k_prev, i_prev2, j_prev2, k_prev2, integrationDirection, dx);
	}

	return P_result;
}

double pressure48(size_t i_dest, size_t j_dest, size_t k_dest, int integrationDirection, int integratonDim, double dx,
float* P_FIELD, float* RHO_FIELD, grid* GRID, bool debug, int level, char* debugLabel)
{
	switch(HydroMethod)
	{
	case Zeus_Hydro:
		return zeusPressure(i_dest, j_dest, k_dest, integrationDirection, integratonDim, dx, P_FIELD, RHO_FIELD, GRID,
							debug, level, debugLabel);
		break;
	case HD_RK:
	case MHD_RK:
	case MHD_Li:
		return rkPressure(i_dest, j_dest, k_dest, integrationDirection, integratonDim, dx, P_FIELD, RHO_FIELD, GRID,
							debug, level, debugLabel);
		break;
	}
	return 0;
}

void PtoE(grid* g, float *totEField, float *pressureField, float* rhofield)
// Note: totEField and pressureField can be the same
{
	const size_t n = g->GetGridSize();
	const double gm1 = Gamma - 1;
	for(size_t i = 0; i < n; i++)
		totEField[i] = pressureField[i] / rhofield[i] / gm1;
}

void set1(grid* g, float* totEField, int i, int j, int k, float P, int i0, int j0, int k0, float* rhoField,
float* vField[3], float* BField[3])
{
	FLOAT LatticeIntegrationRadius = -1;
	FLOAT LatticeIntegrationRadiusInt = 5;
	FLOAT x, y, z, r;

	if(!g->ijkInGrid(i, j, k))
		return;

	if(LatticeIntegrationRadiusInt > 0 || LatticeIntegrationRadius > 0)
	{
		x = g->GetCellCenter(1, i) - SphericalGravityCenter[0];
		y = g->GetCellCenter(1, j) - SphericalGravityCenter[1];
		z = g->GetCellCenter(2, k) - SphericalGravityCenter[2];
		r = lenl(x, y, z);
	}

	if(LatticeIntegrationRadius > 0 && r > LatticeIntegrationRadius)
		return;

	if(LatticeIntegrationRadiusInt > 0)
	{
		r /= g->GetCellWidth(0, 0);
		if(r > LatticeIntegrationRadiusInt)
			return;
	}

	int index = g->GetIndex(i, j, k);
	float rho = rhoField[index];
	float U = square(BField[0][index]) + square(BField[1][index]) + square(BField[2][index]);
	U /= rho;
	U += square(vField[0][index]) + square(vField[1][index]) + square(vField[2][index]);
	U /= 2;
	U += P / rho / (Gamma - 1);
	totEField[index] = U;
	if(i >= 103 && j >= 103 && k >= 103)
		TRACEF("%lld %lld %lld   P=%e  U=%e    rho=%e  Gamma=%e", i, j, k, P, U, rho, Gamma);
}

void set48_6(grid* g, float* totEField, int i, int j, int k, double P, int i0, int j0, int k0, float* rhoField,
float* vField[3], float* BField[3])
{
	set1(g, totEField, i, j, k, P, i0, j0, k0, rhoField, vField, BField);
	set1(g, totEField, i, k, j, P, i0, j0, k0, rhoField, vField, BField);
	set1(g, totEField, j, i, k, P, i0, j0, k0, rhoField, vField, BField);
	set1(g, totEField, j, k, i, P, i0, j0, k0, rhoField, vField, BField);
	set1(g, totEField, k, i, j, P, i0, j0, k0, rhoField, vField, BField);
	set1(g, totEField, k, j, i, P, i0, j0, k0, rhoField, vField, BField);
}

void set48(grid* g, float* totEField, int i, int j, int k, double P, int i0, int j0, int k0, float* rhoField,
float* vField[3], float* BField[3])
{
	int i2 = 2 * i0 - 1 - i;
	int j2 = 2 * j0 - 1 - j;
	int k2 = 2 * k0 - 1 - k;
	set48_6(g, totEField, i, j, k, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i2, j, k, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i, j2, k, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i2, j2, k, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i, j, k2, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i2, j, k2, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i, j2, k2, P, i0, j0, k0, rhoField, vField, BField);
	set48_6(g, totEField, i2, j2, k2, P, i0, j0, k0, rhoField, vField, BField);
}

void iterator48(grid* g, float* field, float* rhoField, MHDInitialProfile* p, int level, float* vField[3],
float* BField[3])
{
	const bool DEBUG_SWEEPS = false;
	const bool DEBUG48 = false;

	const size_t FIELD_SIZE = g->GetGridSize();
	const int I2 = g->GetGridDimension(0);
	const int J2 = g->GetGridDimension(1);
	const int K2 = g->GetGridDimension(2);
	const int I0 = I2 / 2; //TODO: Use SphericalGravityCenter
	const int J0 = J2 / 2;
	const int K0 = K2 / 2;
	const int K1 = K0 + 2;
	const int K0b = K2 - K0 - 1;
	const int K1b = K2 - K1 - 1;
	const FLOAT DX = g->GetCellWidth(0, 0);
	FLOAT x, y, z, r;
	double P_c, P, U;
	double rho, rho_c, ksi, ksifactor;
	int err;
	level = (level >= 0) ? level : nint(log2(TopGridDx[0] / DX));

//	arr_set(field, FIELD_SIZE, 0);

	p->interpolateDensity(&rho_c, 0);

	TRACEF("%lld .... %lld .. %lld | %lld .. %lld .... %lld", 0, K1b, K2 - K0 - 1, K0, K1, K2);

//Compute central zones
	for(size_t k = K0; k < K1; k++)
	{
		for(size_t j = J0; j < K1; j++)
		{
			for(size_t i = I0; i < K1; i++)
			{
				polytropicPressureAtSmallR(&P, NULL, i, j, k, g, rho_c, &P_c, &ksi, &ksifactor, &r);
				P = (DEBUG_SWEEPS) ? 1e1 : P;
				TRACEF("ijk=%lld %lld %lld xyz=%e %e %e xi,r/xi=%e %e   P=%e   P_c=%e  g=%e", i, j, k, x, y, z, ksi,
						ksifactor, P, P_c, SphericalGravityGetAt(r));
				set48(g, field, i, j, k, P, I0, J0, K0, rhoField, vField, BField);
			}
		}
	}

// Integrate from the center cells out along z.
	for(size_t j = J0; j < K1; j++)
	{
		for(size_t i = I0; i < K1; i++)
		{
			for(size_t k = K1; k < K2; k++)
			{
				P = pressure48(i, j, k, 1, 2, DX, field, rhoField, g, (k - i < 6), level, "z/48 sweep");
				P = (DEBUG_SWEEPS) ? 1e2 : P;
				set48(g, field, i, j, k, P, I0, J0, K0, rhoField, vField, BField);
			}
		}
	}
	TRACEF("P_boundary z = %e %e", P, P_c);

// Integrate from the center cells out along y.
	for(size_t k = K1; k < K2; k++)
	{
		for(size_t i = I0; i < K1; i++)
		{
			for(size_t j = K1; j <= k; j++)
			{
				P = pressure48(i, j, k, 1, 1, DX, field, rhoField, g, DEBUG48, level, "y/48 sweep");
				P = (DEBUG_SWEEPS) ? 1e3 : P;
				set48(g, field, i, j, k, P, I0, J0, K0, rhoField, vField, BField);
			}
		}
	}
	TRACEF("P_boundary y = %e %e", P, P_c);

	for(size_t k = 0; k < K2; k++)
	{
		for(size_t j = 0; j <= k; j++)
		{
			for(size_t i = K1; i <= j; i++)
			{
				P = pressure48(i, j, k, 1, 0, DX, field, rhoField, g, DEBUG48, level, "x/48 sweep");
				P = (DEBUG_SWEEPS) ? 1e4 : P;
				set48(g, field, i, j, k, P, I0, J0, K0, rhoField, vField, BField);
			}
		}
	}
	TRACEF("P_boundary x = %e %e", P, P_c);

//	PtoE(g, field, field, rhoField);
}

int grid::MHD_SNIA_GetFields(float** densityField, float** totalEnergyField, float** internalEnergyField,
float** vxField, float** vyField, float** vzField, float** vFields,
float** BxField, float** ByField, float** BzField, float** BFields,
float** burnedDensityField,
float** phiField,
float** gravPotentialField)
{
	int rhoNum, gasENum, vxNum, vyNum, vzNum, totENum;
	int BxNum, ByNum, BzNum, phiNum;
	int rhoNiNum;
	int gravPotentialNum;

// IdentifyPhysicalQuantities set gasENum=0 if not having Gas Energy field.
	bool hasGasEField = DualEnergyFormalism; // && (HydroMethod != Zeus_Hydro); //[BH]
	BxNum = ByNum = BzNum = phiNum = -1; // IdentifyPhysicalQuantities doesn't change these args if (not UseMHD).
	if(IdentifyPhysicalQuantities(rhoNum, gasENum, vxNum, vyNum, vzNum, totENum, BxNum, ByNum, BzNum, phiNum) == FAIL)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in IdentifyPhysicalQuantities for UseMHD==true.")

	if(densityField)
		*densityField = BaryonField[rhoNum];

	if(totalEnergyField)
		*totalEnergyField = BaryonField[totENum];

	if(internalEnergyField)
		*internalEnergyField = (hasGasEField) ? BaryonField[gasENum] : NULL;

	if(vxField)
		*vxField = BaryonField[vxNum];
	if(vyField)
		*vyField = (MaxVelocityIndex > 1) ? BaryonField[vyNum] : NULL;
	if(vzField)
		*vzField = (MaxVelocityIndex > 2) ? BaryonField[vzNum] : NULL;
	if(vFields)
	{
		vFields[0] = BaryonField[vxNum];
		vFields[1] = (MaxVelocityIndex > 1) ? BaryonField[vyNum] : NULL;
		vFields[2] = (MaxVelocityIndex > 2) ? BaryonField[vzNum] : NULL;
	}

	if(UseMHD)
	{
		if(BxField)
			*BxField = BaryonField[BxNum];
		if(ByField)
			*ByField = BaryonField[ByNum];
		if(BzField)
			*BzField = BaryonField[BzNum];

		if(BFields)
		{
			BFields[0] = BaryonField[BxNum];
			BFields[1] = BaryonField[ByNum];
			BFields[2] = BaryonField[BzNum];
		}
	}
	else
	{
		if(BxField)
			*BxField = NULL;
		if(ByField)
			*ByField = NULL;
		if(BzField)
			*BzField = NULL;
		if(BFields)
			BFields[0] = BFields[1] = BFields[2] = NULL;
	}

	if(burnedDensityField)
	{
		if(UseBurning)
		{
			rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
			if(rhoNiNum < 0)
				ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
			*burnedDensityField = BaryonField[rhoNiNum];
		}
		else
		{
			*burnedDensityField = NULL;
		}
	}

	if(phiField)
		*phiField = (phiNum >= 0) ? BaryonField[phiNum] : NULL;

	if(gravPotentialField)
	{
		if(WritePotential)
		{
			gravPotentialNum = FindField(GravPotential, FieldType, NumberOfBaryonFields);
			if(gravPotentialNum < 0)
				ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for GravPotential.")
			*gravPotentialField = BaryonField[gravPotentialNum];
		}
		else
		{
			*gravPotentialField = NULL;
		}
	}

	return SUCCESS;
}

/*
 * Init phase 0
 * Stub the grid without allocating any fields.
 * Used with ParallelRootGridIO.
 */
int grid::MHDProfileInitializeGrid0(MHDInitialProfile* p, float burningTemperature, float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere)
{
	TRACEGF("INITIALIZING GRID PHASE 0 BEGIN.");

	if(GridRank != 3)
		ENZO_FAIL("MHDProfileInitializeGrid is implemented for 3D only.")

	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	if(CellWidth[0][0] <= 0)
		PrepareGridDerivedQuantities();

// Assign fieldType numbers using constants from typedefs.h,
// as well as count the number of fields.
	NumberOfBaryonFields = 0;
	FieldType[NumberOfBaryonFields++] = Density;
	if(EquationOfState == 0)
		FieldType[NumberOfBaryonFields++] = TotalEnergy;
	if(hasGasEField)
		FieldType[NumberOfBaryonFields++] = InternalEnergy;
	FieldType[NumberOfBaryonFields++] = Velocity1;
	FieldType[NumberOfBaryonFields++] = Velocity2;
	FieldType[NumberOfBaryonFields++] = Velocity3;
	if(WritePotential)
		FieldType[NumberOfBaryonFields++] = GravPotential;
	if(UseBurning)
		FieldType[NumberOfBaryonFields++] = Density_56Ni; //[BH]
	if(UseMHD)
	{
		FieldType[NumberOfBaryonFields++] = Bfield1;
		FieldType[NumberOfBaryonFields++] = Bfield2;
		FieldType[NumberOfBaryonFields++] = Bfield3;
		if(HydroMethod == MHD_RK)
			FieldType[NumberOfBaryonFields++] = PhiField;
	}

	bool willUseVectorPotential;
	MHDProfileInitializeGrid_B_and_CT_Fields(NULL, burningTemperature, burnedRadius, dipoleMoment, dipoleCenter,
												useVectorPotential, MetaData, triSphere, NULL, &willUseVectorPotential);

	if(HydroMethod == MHD_RK || willUseVectorPotential)
	{
		// Allow the ElectricField and MagneticField to
		// be created temporarily for the purpose of
		// initializing with vector potential.
		// Free these fields in MHDProfileInitializeGrid5.
		MHD_SetupDims();
	}

	TRACEGF("INITIALIZING GRID PHASE 0 END.");
	return SUCCESS;
}

int grid::PerturbWithTriSPhere(TriSphere *triSphere, FILE *fptr = NULL)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	float *rhoField, *totEField, *vxField, *vyField, *vzField, *vFields[3];
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;
	MHD_SNIA_GetFields(&rhoField, &totEField, NULL, NULL, NULL, NULL, vFields, NULL, NULL, NULL, NULL, &rhoNiField,
	NULL,
						NULL);

	if(fptr)
	{
		fprintf(fptr, "type('', (), dict(\n");

		fprintf(fptr, ""
				"## Unperturbed sphere (primary surface) triangulation parameters:\n"
				"R = %e , # Primary radius\n"
				"nV = %lld , # Number of vertices\n"
				"nE = %lld , # Number of edges\n"
				"nF = %lld , # Number of faucets \n"
				"## Perturbations parameters:\n"
				"A = %e , # Approx. perturbation height\n"
				"bottomBaseSize = %4.2f , # relative to the primary facet (linear)\n"
				"topBaseSize = %4.2f , # relative to the primary facet (linear)\n",
				triSphere->R, triSphere->nV, triSphere->nE, triSphere->nF, triSphere->A, triSphere->bottomBaseSize,
				triSphere->topBaseSize);
		fprintf(fptr, "## Grid parameters:\n"
				"gridLevel = %lld , # Hierarchy level\n"
				"gridID = %lld , # grid ID in level\n"
				"grid_dx = %e , # Grid spatial step\n"
				"numGhostZones = %lld , \n"
				"gridEdges = [\n"
				"  [ %e , %e , %e ], # grid left edges\n"
				"  [ %e , %e , %e ]], #grid right edges\n",
				0, ID, TopGridDx[0], NumberOfGhostZones, GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2],
				GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);

		fprintf(fptr, "\n# Perturbation data for the grid:\ndata=[\n\n");
	}

	//Perturb in the shell between r1 < r <r2.
	// Calculate position with respect to each side.
	for(size_t k_face = 0; k_face < triSphere->nF; k_face++)
	{
		long lijk[3], rijk[3];
		double *ledge, *redge;
		FLOAT xyz[3], xyz2[3];
		size_t numInside = 0;

		triSphere->getSpikeEnclosingRectangle(k_face, &ledge, &redge);
		if(fptr)
		{
			fprintf(fptr, "[ %lld ,\n", k_face);
			triSphere->fprint_cache(fptr);
			fprintf(fptr, "],\n[");
		}

		if(0 > intersect(ledge, redge))
		{
			if(fptr)
				fprintf(fptr, "]],\n");
			continue;
		}

		get_ijk_index(lijk, ledge);
		get_ijk_index(rijk, redge);

		for(size_t k = lijk[2]; k < rijk[2]; k++)
		{
			xyz[2] = CELLCENTER(2, k);
			for(size_t j = lijk[1]; j < rijk[1]; j++)
			{
				xyz[1] = CELLCENTER(1, j);
				for(size_t i = lijk[0]; i < rijk[0]; i++)
				{
					xyz[0] = CELLCENTER(0, i);

					bool isInsideSpike = triSphere->isInSpike(k_face, xyz);
					if(!isInsideSpike)
						continue;

					if(fptr)
						fprintf(fptr, "[ %e , %e , %e ],\n", xyz[0], xyz[1], xyz[2]);

					size_t index = ELT(i, j, k);
					double massfrac = 1;
					//double DX = TopGridDx[0] / CellWidth[0][0];
					//double massfrac = (DX > .75) ? 1 : (DX > .4) ? .2 : .01;
					//massfrac = 4 * massfrac / (1 + 3 * massfrac);

					float rho;
					if(PerturbationBottomDensity > 0)
						rhoField[index] = rho = PerturbationBottomDensity;
					else
					{
						rho = rhoField[index];
						rhoField[index] = rho;
					}

					if(PerturbationVelocity > 0)
					{
						double r = lenl(xyz, 3);
						if(r >= InitialBurnedRadius)
						{
							double dKE = 0;
							for(int dim = 0; dim < GridRank; dim++)
							{
								dKE -= square(vFields[dim][index]);
								vFields[dim][index] = PerturbationVelocity * xyz[dim] / InitialBurnedRadius;
								dKE += square(vFields[dim][index]);
							}
							totEField[index] = dKE / 2;
						}
					}
					rhoNiField[index] = massfrac * rho;
					numInside++;
				}
			}
		} // end i,j,k loops

		if(fptr)
			fprintf(fptr, "]],\n");
	} // end k_face loop

	if(fptr)
	{
		fprintf(fptr, "\n] #end data\n");
		fprintf(fptr, ") #end dict\n");
		fprintf(fptr, ") #end type\n\n");
		fclose(fptr);
		fptr == NULL;
	}

	if(PerturbationVelocity > 0)
	{
		FLOAT xyz[3];
		for(size_t k = 0; k < GridDimension[2]; k++)
		{
			xyz[2] = CELLCENTER(2, k);
			for(size_t j = 0; j < GridDimension[1]; j++)
			{
				xyz[1] = CELLCENTER(1, j);
				for(size_t i = 0; i < GridDimension[0]; i++)
				{
					xyz[0] = CELLCENTER(0, i);
					double r = lenl(xyz, 3);
					if(r >= InitialBurnedRadius)
						continue;

					size_t index = ELT(i, j, k);
					double dKE = 0;
					for(int dim = 0; dim < GridRank; dim++)
					{
						dKE -= square(vFields[dim][index]);
						vFields[dim][index] = PerturbationVelocity * xyz[dim] / InitialBurnedRadius;
						dKE += square(vFields[dim][index]);
					}
					totEField[index] += dKE / 2;
				}
			}
		}
	}

	return SUCCESS;
}

/*
 * Init phase 1
 * Stub (phase 0) if necessary.
 * Allocate fields.
 * Initialize
 *  -- density (eventually from profile);
 *  -- burned fraction;
 * Perturb
 *  -- burned interface;
 *  -- velocity;
 */
int grid::MHDProfileInitializeGrid1(MHDInitialProfile* radialProfile, float burningTemperature, float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(GridRank != 3)
	{
		ENZO_FAIL("MHDProfileInitializeGrid is implemented for 3D only.");
	}

	if(!ParallelRootGridIO)
	{
		if(MHDProfileInitializeGrid0(NULL, burningTemperature, burnedRadius, dipoleMoment, dipoleCenter,
										useVectorPotential, MetaData, triSphere) == FAIL)
		{
			ENZO_FAIL("Error stubbing grid by MHDProfileInitialize0.");
		}
	}

//	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	TRACEGF("INITIALIZING GRID PHASE 1 BEGIN.");

	this->AllocateGrids();

	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField, NULL, &rhoNiField, NULL, NULL);

	int debugnanflag = 0; //[BH]
	float gammaMinusOne = Gamma - 1;
	float rho_c;
	radialProfile->interpolateDensity(&rho_c, 0);

	// In case we want to perturb the front at a later restart,
	// initialize unperturbed burned sphere now (method=0).
	// Methods -3, -2, -1 are unperturbed burned slabs
	// perpendicular to x, y, or z, for test problems.
	// Those are not intendedto be used on restart.
	int perturbMethodNow = (PerturbationOnRestart && (PerturbationMethod > 0)) ? 0 : PerturbationMethod;

// Initialize density and a sphere of Nickel density.
	for(int k = 0; k < GridDimension[2]; k++)
	{
		FLOAT z = CELLCENTER(2, k);
		FLOAT rz = z - SphericalGravityCenter[2];
		FLOAT zz = square(rz);
		for(int j = 0; j < GridDimension[1]; j++)
		{
			FLOAT y = CELLCENTER(1, j);
			FLOAT ry = y - SphericalGravityCenter[1];
			FLOAT yy_zz = square(ry) + zz;
			for(int i = 0; i <= GridDimension[0]; i++)
			{
				int index = ELT(i, j, k);
				float vr, vx, vy, vz, vv, Bx, By, Bz, BB;

				FLOAT x = CELLCENTER(0, i);
				FLOAT rx = x - SphericalGravityCenter[0];
				FLOAT r = sqrt(square(rx) + yy_zz);

				float rho, rhoNi = 0;
				int retcode = radialProfile->interpolateDensity(&rho, r); //g/cm**3
				bool isBurned = false;
				double r1, r2;

				switch(perturbMethodNow)
				{
				case -3: // xy-slab
					isBurned = fabs(rz) < InitialBurnedRadius;
					break;
				case -2: // zx-slab
					isBurned = fabs(ry) < InitialBurnedRadius;
					break;
				case -1: // yz-slab
					isBurned = fabs(rx) < InitialBurnedRadius;
					break;
				case 0:
					isBurned = r <= InitialBurnedRadius;
					break;
				case 1: // Burned sphere, smooth perturb zone.
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
					{
						rhoNi = (r - r1) / PerturbationAmplitude;
						rhoNi *= 4 / (1 + 3 * rhoNi);
						rhoNi *= rho;
					}
					break;
				case 2: // Burned sphere, randomly burned perturb zone
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
						isBurned = rand() % 2;
					break;
				case 3: // Burned sphere, randomly burned perturb zone
						// Probability decreasing with radius.
					r1 = InitialBurnedRadius;
					r2 = r1 + PerturbationAmplitude;
					if(r <= r1)
						isBurned = true;
					else if(r > r2)
						isBurned = false;
					else
					{
						isBurned = drand48() >= (r - r1) / PerturbationAmplitude;
					}
					break;
				case 4:
					isBurned = r <= InitialBurnedRadius;
					break;
				}

				if(isBurned)
				{
//					rho *= 0.75;
					rhoNi = rho;
				}
				rhoField[index] = rho;
				if(rho <= 0)
					TRACEF("Bad density at  %lld %lld %lld :  %e", i, j, k, rho);
				if(UseBurning)
					rhoNiField[index] = rhoNi;

				// Init velocity.
				vr = 0;
				if(radialProfile->radialVelocityData)
					radialProfile->interpolateRadialVelocity(&vr, r);
				else if(0)
				{
					// values for testing
					vr = 10e5; // radial
					if(1) // tangential
					{
						vx = vr * (-ry / r);
						vy = vr * (rx / r);
						vr = 0;
					}
				}
				//Project radial velocity to coordinate components.
				//If the above method sets these components itself
				//it should set vr = 0.
				vx += vr * rx / r;
				vy += vr * ry / r;
				vz += vr * rz / r;

				//Add bulk velocity.
				vx += GridVelocity[0];
				vy += GridVelocity[1];
				vz += GridVelocity[2];
			}
		}
	}

	//Spikes perturbations (method == 4)
	if(perturbMethodNow == 4 && burnedRadius > 0)
	{
		FILE* fptr = NULL;
		if(MyProcessorNumber == ROOT_PROCESSOR && isTopGrid())
			fptr = fopen("trisphere_init.py", "w");
		PerturbWithTriSPhere(triSphere, fptr);
	}

// Boiler plate code:
//  if(DualEnergyFormalism )
//    for(index=0;index<size;index++)
//    BaryonField[ gesENum ][index] =
//       BaryonField[ totENum ][index]
//      - 0.5*(BaryonField[ vNum[0] ][index]*BaryonField[ vNum[0] ][index] +
//	     BaryonField[ vNum[1] ][index]*BaryonField[ vNum[1] ][index] +
//	     BaryonField[ vNum[2] ][index]*BaryonField[ vNum[2] ][index])
//      - 0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index] +
//             BaryonField[ByNum][index]*BaryonField[ByNum][index] +
//             BaryonField[BzNum][index]*BaryonField[BzNum][index])/BaryonField[ rhoNum ][index];

	if(debugnanflag)
		printf("nans intialized.\n"); //[BH]

	TRACEGF("INITIALIZING GRID PHASE 1 END.");
	return SUCCESS;
}

/*
 * Keeps the burned fraction equal to 1 for the initial burned raadius.
 * Using this method is necessary only when initial radius is very small
 * until the front develops.
 */
int grid::MHDSustainInitialBurnedRegionGrid()
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(!UseBurning)
		return SUCCESS;

	int rhoNum = FindField(Density, FieldType, NumberOfBaryonFields);
	if(rhoNum < 0)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
	int rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
	if(rhoNiNum < 0)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")

	float* rhoField = BaryonField[rhoNum];
	float* rhoNiField = BaryonField[rhoNiNum];
	FLOAT lxyz[MAX_DIMENSION], rxyz[MAX_DIMENSION];
	long lindex[MAX_DIMENSION], rindex[MAX_DIMENSION];

	arr_set(lxyz, MAX_DIMENSION, 0);
	arr_set(rxyz, MAX_DIMENSION, 0);
	for(int dim = 0; dim < GridRank; dim++)
	{
		lxyz[dim] = SphericalGravityCenter[dim] - InitialBurnedRadius;
		rxyz[dim] = SphericalGravityCenter[dim] + InitialBurnedRadius;
	}

	if(intersect(lxyz, rxyz))
		return SUCCESS;

	get_ijk_index(lindex, lxyz);
	get_ijk_index(rindex, rxyz);

	const FLOAT RR = square(InitialBurnedRadius);
	const FLOAT X0 = SphericalGravityCenter[0];
	switch(GridRank)
	{
	case 1:
	{
		FLOAT xmin = X0 - RR;
		FLOAT xmax = X0 + RR;
		for(int i = lindex[0]; i <= rindex[0]; i++)
		{
			FLOAT x = CELLCENTER(0, i);
			if(xmin <= x && x <= xmax)
				rhoNiField[i] = rhoField[i];
		}
		break;
	}
	case 2:
		for(int j = lindex[1]; j < rindex[1]; j++)
		{
			FLOAT xxmax = RR - square(CELLCENTER(1, j) - SphericalGravityCenter[1]);
			if(xxmax < 0)
				continue;
			size_t index = ELT(lindex[0], j, 0);
			for(int i = lindex[0]; i <= rindex[0]; i++)
			{
				FLOAT xx = square(CELLCENTER(0, i) - X0);
				if(xx <= xxmax)
					rhoNiField[index] = rhoField[index];
				index++;
			}
		}
		break;

	case 3:
		for(int k = lindex[2]; k <= rindex[2]; k++)
		{
			FLOAT zz = square(CELLCENTER(2, k) - SphericalGravityCenter[2]);
			for(int j = lindex[1]; j < rindex[1]; j++)
			{
				FLOAT xxmax = RR - square(CELLCENTER(1, j) - SphericalGravityCenter[1]) - zz;
				if(xxmax < 0)
					continue;
				size_t index = ELT(lindex[0], j, k);
				for(int i = lindex[0]; i <= rindex[0]; i++)
				{
					FLOAT xx = square(CELLCENTER(0, i) - X0);
					if(xx <= xxmax)
						rhoNiField[index] = rhoField[index];
					// Debug code:
					// rhoNiField[index] = (xx <= xxmax) ? 0 : rhoField[index];
					index++;
				}
			}
		}
		break;
	default:
		ENZO_VFAIL("Bad grid rank: %lld", GridRank)
		;
	}

	return SUCCESS;
}

/*
 * Initialize the gas and total energy from the profile.
 * It might be specified in the input file or calculated by Ezno.
 */
int grid::MHDProfileInitializeGrid_B_and_CT_Fields(MHDInitialProfile* radialProfile, float burningTemperature,
float burnedRadius, float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData,
	TriSphere *triSphere, bool *out_usingDirectInit, bool *out_usingVectorPotentialInit)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	bool direct = true;
	bool vector = false;

	if(out_usingDirectInit)
		*out_usingDirectInit = direct;
	if(out_usingVectorPotentialInit)
		*out_usingVectorPotentialInit = vector;

	if(radialProfile == NULL)
		return SUCCESS;

	TRACEGF("INITIALIZING GRID  E,M,B  BEGIN.");

	MHDClear_B_and_CT_Fields();

	if(direct)
	{
		InitializeMagneticUniformField(BA, 1);
	}

	if(vector)
	{
		InitializeMagneticUniformFieldVectorPotential(BA, 1);
		InitializeMagneticDipoleVectorPotential(dipoleMoment, dipoleCenter, 1);
	}

	TRACEGF("INITIALIZING GRID  E,M,B  END.");
	return SUCCESS;
}

/*
 * It resets the total energy with the updated magnetic field.
 * In the case of using a vector potential fo rintialization or
 * if hydro method  == MHD_RK, it frees the E&M fields.
 * This function is to be used if and after magnetic vector
 * potential has been projected to parents and the curl taken.
 * It is desirable to estimate the magnetic energy in phase 1 or 2
 * because it may affect the intial refinement.
 */
int grid::MHDProfileInitializeGrid_CurlAndCenter(MHDInitialProfile* radialProfile,
float burningTemperature,
float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere,
	bool usingDirectInit, bool usingVectorPotentialInit)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if(!usingDirectInit && !usingVectorPotentialInit)
	{
		MHDProfileInitializeGrid_B_and_CT_Fields(NULL, burningTemperature, burnedRadius, dipoleMoment, dipoleCenter,
													useVectorPotential, MetaData, triSphere, &usingDirectInit,
													&usingVectorPotentialInit);
	}

	if(!usingVectorPotentialInit)
		return SUCCESS;

	TRACEGF("INITIALIZING GRID curl and center BEGIN.");

	int curlMode = (usingDirectInit) ? 3 /*add*/: 0 /*assign*/;
	if(MHD_Curl(curlMode) == FAIL)
		ENZO_FAIL("MHDProfileRefineOnStartup: error in MHD_Curl\n");

	if(CenterMagneticField() == FAIL)
		ENZO_FAIL("MHDProfileRefineOnStartup: error in CenterMagneticField\n");

	TRACEGF("INITIALIZING GRID curl and center END.");
	return SUCCESS;
}

int grid::MHDProfileInitializeGrid2(MHDInitialProfile* radialProfile, float burningTemperature, float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	TRACEGF("INITIALIZING GRID PHASE 2 BEGIN.");

	bool usingDirectInit, usingVectorPotentialInit;

	MHDProfileInitializeGrid_B_and_CT_Fields(radialProfile, burningTemperature, InitialBurnedRadius, dipoleMoment,
//					InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel) ? 1 : 0,
												dipoleCenter, useVectorPotential, MetaData, triSphere, &usingDirectInit,
												&usingVectorPotentialInit);

	MHDProfileInitializeGrid_CurlAndCenter(radialProfile, burningTemperature, InitialBurnedRadius, dipoleMoment,
											dipoleCenter,
//					InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel) ? 1 : 0,
											useVectorPotential, MetaData, triSphere, usingDirectInit,
											usingVectorPotentialInit);

	MHDProfileInitializeGrid_TotalE_GasE(radialProfile, burningTemperature, InitialBurnedRadius, dipoleMoment,
											dipoleCenter,
//					InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel) ? 1 : 0,
											useVectorPotential, MetaData, triSphere);

	TRACEGF("INITIALIZING GRID PHASE 2 END.");
	return SUCCESS;
}

int grid::MHDProfileInitializeGrid_TotalE_GasE(MHDInitialProfile* radialProfile, float burningTemperature,
	float burnedRadius,
	float dipoleMoment[3], float dipoleCenter[3], bool usingVectorPotential, TopGridData *MetaData,
	TriSphere *triSphere)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	TRACEGF("INITIALIZING GRID total/internal energy BEGIN.");

	size_t index;
	float rho, gasE, totE;
	float x, y, z, r, rx, ry, rz;
	float *totEField, *rhoField;
	float *vxField, *vyField, *vzField, *vxyzField[3];
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL; //, *BxyzField[3] = { NULL, NULL, NULL };
	MHD_SNIA_GetFields(&rhoField, &totEField, NULL, &vxField, &vyField, &vzField, NULL, &BxField, &ByField, &BzField,
	NULL,
						NULL, NULL, NULL);
	bool hasInternalEnergyProfile = radialProfile->internalEnergyData;
	const float tiny_gasE = tiny_pressure / (Gamma - 1);

	arr_set(totEField, gridSize, 0);

	for(int k = 0; k < GridDimension[2]; k++)
	{
		z = CELLCENTER(2, k);
		rz = z - SphericalGravityCenter[2];

		for(int j = 0; j < GridDimension[1]; j++)
		{
			y = CELLCENTER(1, j);
			ry = y - SphericalGravityCenter[1];
			index = ELT(0, j, k);

			for(int i = 0; i <= GridDimension[0]; i++)
			{
				x = CELLCENTER(0, i);
				rx = x - SphericalGravityCenter[0];
				r = lenl(rx, ry, rz);

				rho = max(rhoField[index], tiny_number);

				if(hasInternalEnergyProfile)
					radialProfile->interpolateInternalEnergy(&gasE, r);
				else
					gasE = internalEnergy(rho, rhoNiField[index], radialProfile, r);

				gasE = max(gasE, tiny_gasE);

				totE = 0;
				if(BxField)
				{
					totE += square(BxField[index]) + square(ByField[index]) + square(BzField[index]);
					totE /= rho;
				}
				totE += square(vxField[index]) + square(vyField[index]) + square(vzField[index]);
				totE /= 2;
				totE += gasE;
				totEField[index] = totE;
				if(gasEField)
					gasEField[index] = gasE;

				index++;
			}
		}
	}

	// Boiler plate code:
	//  if(DualEnergyFormalism )
	//    for(index=0;index<size;index++)
	//    BaryonField[ gesENum ][index] =
	//       BaryonField[ totENum ][index]
	//      - 0.5*(BaryonField[ vNum[0] ][index]*BaryonField[ vNum[0] ][index] +
	//	     BaryonField[ vNum[1] ][index]*BaryonField[ vNum[1] ][index] +
	//	     BaryonField[ vNum[2] ][index]*BaryonField[ vNum[2] ][index])
	//      - 0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index] +
	//             BaryonField[ByNum][index]*BaryonField[ByNum][index] +
	//             BaryonField[BzNum][index]*BaryonField[BzNum][index])/BaryonField[ rhoNum ][index];

	TRACEGF("INITIALIZING GRID total/internal energy END.");
	return SUCCESS;
}

int grid::MHDProfileInitializeGrid5(MHDInitialProfile* radialProfile, float burningTemperature, float burnedRadius,
float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere)
{
	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	TRACEGF("INITIALIZING GRID  clean up  BEGIN.");

	if(HydroMethod == MHD_RK)
	{
		for(int dim = 0; dim < 3; dim++)
		{
			delete MagneticField[dim];
			delete ElectricField[dim];
			MagneticField[dim] = NULL;
			ElectricField[dim] = NULL;
		}
	}

	TRACEGF("INITIALIZING GRID  clean up  END.");
	return SUCCESS;
}

//int grid::MHDProfileInitializeGrid6(MHDInitialProfile* radialProfile, float burningTemperature, float burnedRadius,
//float dipoleMoment[3], float dipoleCenter[3], bool useVectorPotential, TopGridData *MetaData, TriSphere *triSphere)
//{
//	if(ProcessorNumber != MyProcessorNumber)
//		return SUCCESS;
//
//	TRACEGF("INITIALIZING GRID PHASE 6 BEGIN.");
//	float *totEField, *rhoField;
//	float *vxField, *vyField, *vzField, *vxyzField[3];
//	float *gasEField = NULL, *rhoNiField = NULL;
//	float *BxField = NULL, *ByField = NULL, *BzField = NULL; //, *BxyzField[3] = { NULL, NULL, NULL };
//	MHD_SNIA_GetFields(&rhoField, &totEField, NULL, NULL, NULL, NULL, NULL, &BxField, &ByField, &BzField, NULL, NULL,
//	NULL, NULL);
//
//	TRACEGF("INITIALIZING GRID PHASE 6 END.");
//	return SUCCESS;
//}

int grid::WriteRadialProfile(char* name)
{
	if(MyProcessorNumber != ProcessorNumber)
		return FAIL;

	long num, ijk[3];
	FLOAT x, y, z, r;
	float rho, E, U, vx, vy, vz;
	float Bx = 0, By = 0, Bz = 0, rho_Ni = 0, g = 0, gamma = 0, P = 0;
	char filename[FILENAME_MAX];
	FLOAT xyz1[MAX_DIMENSION], xyz2[MAX_DIMENSION];
	FLOAT dx = CellWidth[0][0];
	size_t i1 = 0, i2 = 0;

	const int SAMPLE_OX = 0;
	const int SAMPLE_DIAG = 1;
	int sampleWhere = SAMPLE_DIAG;
	switch(sampleWhere)
	{
	case SAMPLE_OX:
// Sample a row of cells along x,
// (0,0,0)..(DomainRightEdge[0], 0, 0)
		arr_set(xyz1, MAX_DIMENSION, dx * 1e-6);
		arr_set(xyz2, MAX_DIMENSION, dx * 1e-6);
		xyz2[0] = DomainRightEdge[0] - dx * 1e-6;
		break;
	case SAMPLE_DIAG:
	default:
// Sample cells from the center along the diagonal,
// (0,0,0)..(DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2]):
		arr_set(xyz1, MAX_DIMENSION, dx * 1e-6);
		arr_set(xyz2, MAX_DIMENSION, -dx * 1e-6);
		arr_xpy(xyz2, DomainRightEdge, MAX_DIMENSION);
	}

//Check if the segment xyz1..xyz2 is inside the grid.
	for(int dim = 0; dim < GridRank; dim++)
	{
		if(xyz1[dim] < GridLeftEdge[dim] || (GridRightEdge[dim] < xyz2[dim]))
			return FAIL;
	}

	switch(sampleWhere)
	{
	case SAMPLE_OX:
	case SAMPLE_DIAG:
	default:
		get_ijk_index(ijk, xyz1);
		i1 = ijk[0];
		get_ijk_index(ijk, xyz2);
		i2 = ijk[0];
	}

	int hasGasEField = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField,
						NULL,
						&rhoNiField, NULL, NULL);

	for(int format = 0; format < 2; format++)
	{
		switch(format)
		{
		case 0:
			sprintf(filename, "%s.radialProfile", name);
			break;
		case 1:
			sprintf(filename, "%s.radialProfile.py", name);
			break;
		}

		FILE* file = fopen(filename, "w");
		if(file == NULL)
		{
			ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
		}

		switch(format)
		{
		case 0:
			fprintf(file, "%f seconds\n", Time);
			break;
		case 1:
			break;
		}

		switch(format)
		{
		case 0:
			fprintf(file, "# rownum   x y z   rho rho_Ni   E U   vx vy vz   Bx By Bz   g   gamma P\n");
			break;
		case 1:
			fprintf(file, "type('', (), dict(\n");
			fprintf(file, "time = %f,\n", Time);
			fprintf(file, "cols = 'rownum x y z rho rho_Ni E U vx vy vz Bx By Bz g gamma P'.split(),\n");
			fprintf(file, "data = np.array([\n");
			break;
		}

		arr_set(ijk, MAX_DIMENSION, 0);
		for(size_t i = i1; i <= i2; i++)
		{
			num = i - i1 + 1;

			switch(sampleWhere)
			{
			case SAMPLE_OX:
				ijk[0] = i;
				break;
			case SAMPLE_DIAG:
			default:
				arr_set(ijk, GridRank, i);
			}

			size_t index = getCellIndex(ijk);

			get_xyz(xyz1, ijk);
			x = xyz1[0];
			y = xyz1[1];
			z = xyz1[2];
			r = lenl(xyz1, SphericalGravityCenter, GridRank);
			E = totEField[index];
			rho = rhoField[index];
			vx = vxField[index];
			vy = vyField[index];
			vz = vzField[index];
			gamma = Gamma;

			if(BxField)
			{
				Bx = BxField[index];
				By = ByField[index];
				Bz = BzField[index];
			}

			if(gasEField)
				U = gasEField[index];
			else
			{
				U = square(Bx) + square(By) + square(Bz);
				U /= rho;
				U += square(vx) + square(vy) + square(vz);
				U *= -0.5;
				U += E;
			}
			P = (gamma - 1) * U * rho;

			if(UseBurning)
				rho_Ni = rhoNiField[index];

			if(UseSphericalGravity)
			{
				g = SphericalGravityGetAt(r);
			}

//		fprintf(file, "num   x   y   z   rho   E   U   vx   vy   vz   Bx   By   Bz   g   gamma P\n");
			switch(format)
			{
			case 0:
				fprintf(file, "%lld   %e %e %e   %e %e   %e %e   %e %e %e   %e %e %e   %e   %f %e\n", /**/
						num, x, y, z, rho, rho_Ni, E, U, vx, vy, vz, Bx, By, Bz, g, gamma, P);
				break;
			case 1:
				fprintf(file, "[ %lld ,   %e , %e , %e ,   %e , %e ,   %e , %e ,   %e , %e , %e ,   %e , %e , %e"
						" ,   %e ,   %f , %e ],\n",
						num, x, y, z, rho, rho_Ni, E, U, vx, vy, vz, Bx, By, Bz, g, gamma, P);
				break;
			}
		}

		switch(format)
		{
		case 0:
			break;
		case 1:
			fprintf(file, "])))\n");
			break;
		} // end for data

		if(file != stderr && file != stdout)
		{
			fclose(file);
			fprintf(stderr, "'%s' written.\n", filename);
		}
	} // end for format

	return SUCCESS;
}
