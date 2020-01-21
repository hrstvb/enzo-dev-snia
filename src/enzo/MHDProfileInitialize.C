/***********************************************************************
 *
 *  Initializes the MHD problem (ProblemType 501) from external radial
 *  profiles.
 *
 *  written by: boyan hristov
 *              based on MHDBlastInititalize by David Collins
 *  date:       2018-
 *  modified1:
 *
 *  PURPOSE:  Problem initializer
 *
 *  PARAMETERS:
 *
 *	ProblemType = 501
 *       ProfileFileName = "myprofiledata.txt"
 *       ProfileFormat = "PLAIN" or "PAH1" or "PAH2"
 *       ProfileType = "RADIAL"
 *       RadiusColumnName = "radius"
 *       DensityColumnName = "density"
 *       TemperatureColumnName = "temperature"
 *       BurnedRadius = 1e5
 *       BurningTemperature = 8.5e9
 *
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "DebugMacros.h"
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelArrayIterator.h"
#include "DebugTools.h"
#include "TriSpherePerturb.C"

class ExternalBoundary;

/* This initializes genearlized discontinuities.  Good for blast waves, Kelving Helmholtz, shock tubes.  Probably others.

 Code flow:
 1.) Declare/Define Defaults for  parameters.
 2.) Read parameters from file.
 3.) Calculate TotalEnergy from quantity given (options are Total Energy, Gas Energy, Pressure.)
 4.) Set up data labels, units.
 5.) Declare Hierarchy Object.
 6.) Define linked list
 7.) Call initializer on each level.
 8.) Project to parent.

 */

void limiter1nInit();
int MHDCTSetupFieldLabels();
int SphericalGravityComputePotential(LevelHierarchyEntry *LevelArray[], TopGridData* MetaData,
	bool interpolateMissingShells);
float SphericalGravityGetAt(FLOAT r);
int polytropicPressureAtSmallR(float* P, float* rho, FLOAT r, float rho_c, float * P_c = NULL, double* ksi = NULL,
	double* ksiFactor = NULL);

void WriteInitialProfile(const char* name, FLOAT* RR, FLOAT* RHO, FLOAT* GG, FLOAT* PP, FLOAT* UU, size_t n,
FLOAT gamma)
{
	if(MyProcessorNumber != ROOT_PROCESSOR)
		return;

	const int M = 10;
	size_t num;
	FLOAT r, dr;
	float rho, P, U, g;
	char filename[FILENAME_MAX];
	sprintf(filename, "%s.initialProfile", name);
	FILE* file = fopen(filename, "w");
	if(file == NULL)
	{
		ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
	}
	fprintf(file,
			"# rownum r                dr               rho              P                U                g            v_r\n");
	for(size_t i = 0; i < n; i++)
	{
		num = i + 1;
		r = RR[i];
		dr = (i) ? (r - RR[i - 1]) : 0;
		rho = RHO[i];
		P = PP[i];
		U = UU[i];
		g = GG[i];
		fprintf(file, "%6lld   %12.9e  %12.9e  %12.9e  %12.9e  %12.9e  %12.9e\n", num, r, dr, rho, P, U, g);
	}
	if(file != stderr && file != stdout)
	{
		fclose(file);
		fprintf(stderr, "'%s' written.\n", filename);
	}

	sprintf(filename, "%s.initialProfile.py", name);
	file = fopen(filename, "w");
	if(file == NULL)
	{
		ENZO_VFAIL("WriteRadialProfile: Error creating file '%s'\n", filename);
	}

	fprintf(file, "type('', (), dict(\n");
	fprintf(file, "time = %f,\n", 0);
	fprintf(file, "gamma = %f,\n", gamma);
//	fprintf(file, "K = %12.9e,\n", K);
	fprintf(file, "cols = 'r dr rho P U g'.split(),\n");

	fprintf(file, "%s = np.array([", "r");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", RR[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "dr");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", (i) ? (RR[i] - RR[i - 1]) : 0);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "rho");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", RHO[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "P");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", PP[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "U");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", UU[i]);
	}
	fprintf(file, "\n]),\n");

	fprintf(file, "%s = np.array([", "g");
	for(size_t i = 0; i < n; i++)
	{
		if(i % M == 0)
			fprintf(file, "\n");
		fprintf(file, "%12.9e,  ", GG[i]);
	}
	fprintf(file, "\n])\n))\n");

	if(file != stderr && file != stdout)
	{
		fclose(file);
		fprintf(stderr, "'%s' written.\n", filename);
	}

}

int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z)
{
	*Bx = BA[0];
	*By = BA[1];
	*Bz = BA[2];
	return SUCCESS;
}

int MHDProfileInitExactB(float B[3], FLOAT x[3])
{
	return MHDProfileInitExactB(B, B + 1, B + 2, x[0], x[1], x[2]);
}

TriSphere* NewTriSphereFromParams(int maxRefLevel)
{
	//		if(0)
	//		{
	//			//Test instances of the Tri sphere
	//			triSphere = new TriSphere(InitialBurnedRadius, PerturbationAmplitude, -0, 1.0, 0.0);
	//			triSphere->writePythonFile("triSphere-data-0-10-00.py");
	//			triSphere->~TriSphere();
	//			ENZO_FAIL("FINISHED TRISPHERE TEST OUTPUT");
	//		}

	TriSphere *ts = NULL;
	if(PerturbationWavelength > 0)
	{
		// A positive PerturbationWaveLength leads to a trisphere with facets with
		// approx. size = dx * PerturbationWavelength,
		// where dx is the resolution of the finest refinement level.
		FLOAT dx = TopGridDx[0] / pow(2, maxRefLevel);
		ts = new TriSphere(InitialBurnedRadius, PerturbationAmplitude, dx, PerturbationWavelength,
							PerturbationBottomSize, PerturbationTopSize);
	}
	else
	{
		// A negative PerturbationWaveLength leads to a trisphere with
		// 20 * 4**(PerturbationWaveLength) facets.
		ts = new TriSphere(InitialBurnedRadius, PerturbationAmplitude, (size_t) (-PerturbationWavelength),
							PerturbationBottomSize, PerturbationTopSize);
	}
	if(MyProcessorNumber == ROOT_PROCESSOR)
		ts->writePythonFile("trisphere.py");

	return ts;
}

TriSphere* NewTriSphereFromParams()
{
	int maxRefLevel = max(0, min(MaximumRefinementLevel, MAX_DEPTH_OF_HIERARCHY));
	return NewTriSphereFromParams(maxRefLevel);
}

int MHDProfileSetLabels(bool useGE, bool UseMetal)
{
	char *DensName = "Density";
	char *TEName = "TotalEnergy";
	char *Vel1Name = "x-velocity";
	char *Vel2Name = "y-velocity";
	char *Vel3Name = "z-velocity";
	char *GPotName = "GravPotential";
	char *BxName = "Bx";
	char *ByName = "By";
	char *BzName = "Bz";
	char *PhiName = "Phi";
	char *GEName = "GasEnergy";
	char *MetalName = "Metal_Density";
	char *Density_56NiName = "Density_56Ni";   //[BH]
	char *MHDCT_BxName = "BxF";
	char *MHDCT_ByName = "ByF";
	char *MHDCT_BzName = "BzF";
	char *MHDCT_ExName = "Ex";
	char *MHDCT_EyName = "Ey";
	char *MHDCT_EzName = "Ez";

	// Ignore the unit strings.
	arr_set(DataUnits, MAX_NUMBER_OF_BARYON_FIELDS, NULL);

	int i = 0;
	DataLabel[i++] = DensName; //"g/cm**3";

	if(EquationOfState == 0)
		DataLabel[i++] = TEName; //"erg/g";

	if(useGE)
		DataLabel[i++] = GEName; //"erg/g";

	DataLabel[i++] = Vel1Name; //"cm/s";
	DataLabel[i++] = Vel2Name; //"cm/s";
	DataLabel[i++] = Vel3Name; //"cm/s";

	if(UseMetal)
	{
		DataLabel[i++] = MetalName;
		NColor++;
	}

	if(WritePotential)
		DataLabel[i++] = GPotName;

	if(UseBurning)
	{
		DataLabel[i++] = Density_56NiName; //[BH]
		NColor++;
	}

	if(UseMHD)
	{
		DataLabel[i++] = BxName;
		DataLabel[i++] = ByName;
		DataLabel[i++] = BzName;
		if(HydroMethod == MHD_RK)
			DataLabel[i++] = PhiName;
	}
	if(UseMHDCT)
		MHDCTSetupFieldLabels();

	return SUCCESS;
}

int SphericalGravityIntegratePressure(MHDInitialProfile* profile, float initialPressure, bool initialPressureAtSurface)
{
	bool fromProfile = !UseSphericalGravity;

	size_t N = (fromProfile) ? profile->nRows : SphericalGravityActualNumberOfBins;
	FLOAT* rr = (fromProfile) ? profile->radiusData : SphericalGravityBinLeftEdges;
	FLOAT* PP = new float[N];
	FLOAT* UU = new float[N];
	FLOAT* mm = (fromProfile) ? (new FLOAT[N]) : NULL;
	FLOAT* dd = new FLOAT[N];
	FLOAT* gg = new FLOAT[N];
	FLOAT* ggdd = new FLOAT[N];

	FLOAT P_c, P2, P1, P, r2, r1, r, U2, U1, U, m2, m1, m; // r2 < r1 < r
	FLOAT h, b, a, bb, aa, g, rho, g_rho, rho_c; // b:=r-r2 > a:=r-r1
	int retcode;
	int j;

	retcode = profile->interpolateDensity(&rho_c, 0); //g/cm**3
	polytropicPressureAtSmallR(NULL, NULL, rr[0], rho_c, &P_c, NULL, NULL);
	for(j = 0; j < N; j++)
	{
		retcode = profile->interpolateDensity(&rho, rr[j]); //g/cm**3
		dd[j] = rho;
	}

	// Integrate the mass
//	m1 = r1 = 0;
//	if(gg)
//	{
//		for(size_t i = 0; i < N; i++)
//		{
//			r = rr[i];
//			if(r == 0)
//			{
//				m = g = 0;
//			}
//			else if(fromProfile)
//			{
//				h = r - r1;
//				profile->interpolateDensity(&rho, r);
//				m = m1 + h * 4 * M_PI * r * r * rho;
//				g = SphericalGravityConstant * m / r / r;
//			}
//			else
//			{
//				g = SphericalGravityGetAt(r);
//			}
//			if(mm)
//				mm[i] = m;
//			if(gg)
//				gg[i] = g;
//
//			r1 = r;
//			m1 = m;
//		}
//	}
	for(j = 0; j < N; j++)
		gg[j] = -SphericalGravityGetAt(rr[j]);

	for(j = 0; j < N; j++)
		ggdd[j] = (gg[j] + 0 * rr[j] + 0) * dd[j];

	for(j = 0; j < N; j++)
		PP[j] = 0;

	int integrationStartIndex = 0;

	if(0)
	{
		int k, minj;
		double miny;
		double a, b, a2, b2, c1, c2, z, z2;
		double zz[N], zz2[N];
		for(j = 2; j < N - 2; j++)
		{
			c1 = 1;
			c2 = z2 = 0;
			r = rr[j];

			k = 1;
			a = r - rr[j - k];
			b = rr[j + k] - r;
			z = b * ggdd[j - k] + a * ggdd[j + k] - (a + b) * ggdd[j];
			z /= a * b * (a + b) / 2;
			zz[j] = z;

			c1 = a - b;
			if(fabs(c1) > 1e-8)
			{

				k = 2;
				a = r - rr[j - k];
				b = rr[j + k] - r;
				c2 = a - b;
				if(fabs(c2) > 1e-8 && fabs(c1 - c2) > 1e-8)
				{
					z = b * ggdd[j - k] + a * ggdd[j + k] - (a + b) * ggdd[j];
					z /= a * b * (a + b) / 2;
					zz2[j] = (c1 * zz[j] - c2 * z) / (c1 - c2);
				}
				else
				{
					c1 = c2 = 0;
					zz2[j] = zz[j];
				}
			}
			else
			{
				c1 = c2 = 0;
				zz2[j] = zz[j];
			}
		}
	}

//	TRACEF("  integrationStartIndex = %lld", integrationStartIndex);
	if(integrationStartIndex)
	{
		// Backward integration from start index to 0.
		j = integrationStartIndex;
		r2 = rr[j];
		PP[j] = P2 = 0;

		j--;
		r1 = rr[j];
		g_rho = ggdd[j];
		a = r1 - r2;
		PP[j] += P1 = P2 + a * g_rho;
		TRACEF("  %lld -> 0", j);
		for(j--; j >= 0; j--)
		{
			r = rr[j];
			b = r - r2;

			a = r - r1;
			bb = square(b);
			aa = square(a);

			g_rho = ggdd[j];
			P = (bb * (P1 + a * g_rho) - aa * (P2 + b * g_rho)) / (bb - aa);
//			P = P1 + a * g * rho;
			// TRACEF("%lld: P2 P1 P = %e %e %e ; r2 r1 r = %e %e %e ; g rho = %e %e", i, P2, P1, P, r2, r1, r, g, rho);

			PP[j] += P;

			r2 = r1;
			r1 = r;
			P2 = P1;
			P1 = P;
		}

		j = integrationStartIndex;
		r2 = rr[j - 1];
		P2 = PP[j - 1];
		r1 = rr[j];
		P1 = PP[j];
	}
	else
	{
		j = integrationStartIndex;
		// Start from the certer with P_c = 0, apply the initial condition later.
		P_c = 0;
		j = 0;
		r = r2 = rr[j];
		polytropicPressureAtSmallR(&P2, NULL, r2, rho_c, &P_c, NULL, NULL);
		PP[j] += P2;
		//TRACEF("r[0], P[0], g[0]rho[0], g[0], rho[0] = %e, %e, %e = %e * %e", r2, P2, g_rho, g, rho);

		j = 1;
		r = r1 = rr[j];
		polytropicPressureAtSmallR(&P1, NULL, r1, rho_c, &P_c, NULL, NULL);
		PP[j] += P1;
		// TRACEF("r[1], P[1], dP[1], g, rho = %e, %e, %e = %e * %e", r1, P1, P1 - P2, g, rho);
	}

	// Forward integration from forward index to N.
//	TRACEF("  %lld -> %lld", j, N);
	for(j++; j < N; j++)
	{
		r = rr[j];
		b = r - r2;
		a = r - r1;
		bb = square(b);
		aa = square(a);

		g_rho = ggdd[j];
		P = (bb * (P1 + a * g_rho) - aa * (P2 + b * g_rho)) / (bb - aa);
//		P = P1 + a * g * rho;
		// TRACEF("%lld: P2 P1 P = %e %e %e ; r2 r1 r = %e %e %e ; g rho = %e %e", i, P2, P1, P, r2, r1, r, g, rho);

		PP[j] += P;

		r2 = r1;
		r1 = r;
		P2 = P1;
		P1 = P;
	}

	// Add the boundary condition (P_shift below) to PP[].
	double P_shift = initialPressure;
	if(P_shift <= 0)
		polytropicPressureAtSmallR(&P_shift, NULL, rr[0], rho_c, NULL, NULL, NULL);
	;
	size_t initialConditionIndex = (initialPressureAtSurface) ? (N - 1) : 0;
//	TRACEF("P_shift=%e %e %e   %lld %lld %e", P_shift, initialPressure, P_c, initialConditionIndex,
//			initialPressureAtSurface, PP[initialConditionIndex]);
	P_shift -= PP[initialConditionIndex];

	// TRACEF("Pressure shift: %e = %e - %e", P2, P1, P);
	// Shift pressure and compute the gas energy.
	for(size_t i = 0; i < N; i++)
	{
		P = PP[i] + P_shift;
		r = rr[i];
		rho = dd[i]; //g/cm**3
		U = P / rho / (Gamma - 1);
		PP[i] = P;
		UU[i] = U;
	}

	profile->internalEnergyData = UU;
	profile->internalEnergyRadiusData = rr;
	profile->nRowsInternalEnergy = N;

	if(MyProcessorNumber == ROOT_PROCESSOR)
		WriteInitialProfile("initialProfile", rr, dd, gg, PP, UU, N, /*K,*/Gamma);

	if(MyProcessorNumber == ROOT_PROCESSOR)
	{
		TRACEF("  Number of points after pressure integration: Profile: %lld,  Internal energy: %lld,  Spherical gravity: %lld",
			   profile->nRows, N, SphericalGravityActualNumberOfBins);
		int n = 100;
		n = (n < N) ? N : n;
//		arr_fprintpylist(stderr, profile->radiusData, n, "radius=");
//		arr_fprintpylist(stderr, profile->densityData, 100, "density=");
//		arr_fprintpylist(stderr, profile->internalEnergyData, 100, "internalEnergy=");
//		arr_fprintpylist(stderr, PP, n, "pressure=");
	}

	delete gg;
	delete mm;
	delete PP;
	//rr and UU are referenced by profile->;

	return SUCCESS;
}

int MHDProfileInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
	ExternalBoundary &Exterior, bool BeforeGridDistribution)
{
	if(MyProcessorNumber != ROOT_PROCESSOR)
		return SUCCESS;

	if(ParallelRootGridIO && BeforeGridDistribution)
	{
		// Stub the top grid, i.e. do not allocate the fields.
		return TopGrid.GridData->MHDProfileInitializeGrid1(NULL, 0, 0, NULL, NULL, 0, &MetaData, NULL);
	}

	// Parameters and their defaults.
	char line[MAX_LINE_LENGTH];
	int UseMetal = FALSE;
	int useGE = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]

	char ProfileFileName[MAX_LINE_LENGTH];
	char ProfileFormat[16]; //"PLAIN" or "PAH1" or "PAH2"
	char ProfileType[16]; //"RADIAL"
	char RadiusColumnName[16]; //"radius"
	char RadialVelocityColumnName[16]; //"v_r"
	char DensityColumnName[16]; //"density"
	char InternalEnergyColumnName[16];
	char TemperatureColumnName[16]; //"temperature"
	float BurningTemperature = 0;
	float ProfileAtTime = -1;
	int ProfileUseFrameTime = 0;
	float dipoleMoment[3] = { 0, 0, 0 };
	float dipoleCenter[3] = { 0, 0, 0 };

	// For MHD_RK and HD_RK
	NSpecies = NColor = 0;
	NoMultiSpeciesButColors = 1;

	if(MHD_LI_GRAVITY_AFTER_PLMPRED)
	{
		fprintf(stderr, "Using modified MHD_Li solver (see MHD_LI_GRAVITY_AFTER_PLMPRED)\n");
	}

	*ProfileFileName = *ProfileFormat = *ProfileType = '\0';
	*RadiusColumnName = *DensityColumnName = *RadialVelocityColumnName = *TemperatureColumnName = '\0';
	*InternalEnergyColumnName = '\0';

	int ret = 0, ret2 = 0;
	while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
	{
		ret = 0;
//		Use PSYM or FSYM for floats, ISYM for ints.
		ret += sscanf(line, "ProfileFileName = %s", &ProfileFileName);
		ret += sscanf(line, "ProfileType = %15s", &ProfileType);
		ret += sscanf(line, "ProfileFormat = %15s", &ProfileFormat);
		ret += sscanf(line, "RadiusColumnName = %15s", &RadiusColumnName);
		ret += sscanf(line, "RadialVelocityColumnName = %15s", &RadialVelocityColumnName);
		ret += sscanf(line, "DensityColumnName = %15s", &DensityColumnName);
		ret += sscanf(line, "InternalEnergyColumnName = %15s", &InternalEnergyColumnName);
		ret += sscanf(line, "TemperatureColumnName = %15s", &TemperatureColumnName);
		ret += sscanf(line, "BurningTemperature = %"FSYM, &BurningTemperature);
		ret += sscanf(line, "ProfileAtTime = %"FSYM, &ProfileAtTime);
		ret += sscanf(line, "ProfileUseFrameTime = %"ISYM, &ProfileUseFrameTime);

		ret += sscanf(line, "BA = %"FSYM" %"FSYM" %"FSYM, BA, BA + 1, BA + 2);

		ret += sscanf(line, "MagneticDipoleMoment = %"FSYM" %"FSYM" %"FSYM, dipoleMoment, dipoleMoment + 1,
						dipoleMoment + 2);
		ret += sscanf(line, "MagneticDipoleCenter = %"FSYM" %"FSYM" %"FSYM, dipoleCenter, dipoleCenter + 1,
						dipoleCenter + 2);
	}

	//M_2_SQRTPI = 2/sqrt(pi)
	const double ONE_OVER_SQRT_4PI = M_2_SQRTPI / 4.0; // = 1/sqrt(2*pi)
	arr_ax(BA, 3, ONE_OVER_SQRT_4PI);
	arr_ax(dipoleMoment, 3, ONE_OVER_SQRT_4PI);

	MHDProfileSetLabels(useGE, UseMetal);
//	if(ParallelRootGridIO && BeforeGridDistribution)
//	{
//		// Stub the top grid, i.e. do not allocate the fields.
//		return TopGrid.GridData->MHDProfileInitializeGrid1(NULL, BurningTemperature, InitialBurnedRadius, dipoleMoment,
//															dipoleCenter, InitBWithVectorPotential, &MetaData);
//	}
//
	if(BurningDiffusionRateReduced <= 0)
		BurningDiffusionRateReduced = BurningDiffusionRate / (DomainRightEdge[0] - DomainLeftEdge[0])
				* MetaData.TopGridDims[0];
	if(BurningReactionRateReduced <= 0)
		BurningReactionRateReduced = BurningReactionRate * (DomainRightEdge[0] - DomainLeftEdge[0])
				/ MetaData.TopGridDims[0];

	if(SphericalGravityOuterRadius <= 0)
	{
		// If outer radius is not specified, set it to the minimal
		// radius to cover the entire domain, i.e. the distance
		// from the gravity center to the farthest domain vertice.
		SphericalGravityOuterRadius = distMaxl(DomainLeftEdge, DomainRightEdge, SphericalGravityCenter,
												MetaData.TopGridRank);
	}

	MHDInitialProfile radialProfile = MHDInitialProfile();
	radialProfile.init(RadiusColumnName, DensityColumnName, InternalEnergyColumnName, TemperatureColumnName,
						RadialVelocityColumnName, DensityProfileMaxRadius, DensityProfileMinDensity);
	radialProfile.read(ProfileFileName, ProfileFormat, ProfileAtTime);

	if(ProfileUseFrameTime)
	{
		TopGrid.GridData->SetTime(MetaData.Time = radialProfile.frameTimeFound);
		printf("Initial time from profile frame: %f", radialProfile.frameTimeFound);
	}

	TriSphere *triSphere = (PerturbationMethod == 4) ? NewTriSphereFromParams() : NULL;
//		arr_cpy(StaticRefineShellCenter[0], SphericalGravityCenter, MAX_DIMENSION);
//		StaticRefineShellInnerRadius[0] = InitialBurnedRadius + min(0, PerturbationAmplitude);
//		StaticRefineShellOuterRadius[0] = InitialBurnedRadius + max(PerturbationAmplitude, 0);
//		StaticRefineShellLevel[0] = MAX_DEPTH_OF_HIERARCHY;
//		StaticRefineShellWithBuffer[0] = 0;

	int maxRefLevel = max(0, min(MaximumRefinementLevel, MAX_DEPTH_OF_HIERARCHY));
	maxRefLevel *= (bool) RefineOnStartup;
	int maxGravityLevel = max(0, min(SphericalGravityMaxHierarchyLevel, maxRefLevel));
	bool projectChildrenToParents = maxRefLevel;
	bool initBWithVectorPotential = InitBWithVectorPotential || (dipoleMoment[0] || dipoleMoment[1] || dipoleMoment[2]);
	bool integratePressure = (radialProfile.internalEnergyData == NULL) || InitRadialPressureFromCentral;

	RebuildHierarchyIterator rhit = RebuildHierarchyIterator(maxRefLevel, &TopGrid, &MetaData);
	for(grid* g = rhit.first(); g; g = rhit.next())
	{
		if(rhit.startingNewLevel)
			TRACEF("Refining hierarchy, level %lld", rhit.currentLevel);

		// Initialize the density and the burned fraction, velocity and perturnb.
		g->MHDProfileInitializeGrid1(&radialProfile, BurningTemperature,
										InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel),
//										InitialBurnedRadius,
										dipoleMoment, dipoleCenter, initBWithVectorPotential,
										&MetaData, triSphere);

		if(!rhit.isLastOnCurrentLevel())
			continue;

		PRINT_HIERARCHY0(rhit.levelArray);

		if(UseSphericalGravity)
		{
			// Update the gravity potential with the new layer density.
			SphericalGravityComputePotential(rhit.levelArray, &MetaData, true);

			if(integratePressure)
				SphericalGravityIntegratePressure(&radialProfile, InitRadialPressureFromCentral, false);
		}

		LevelArrayIterator it2 = LevelArrayIterator(rhit.levelArray);
		for(grid* g2 = it2.firstFromLevel(rhit.currentLevel); g2; g2 = it2.next())
			g2->MHDProfileInitializeGrid2(&radialProfile, BurningTemperature,
											InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel),
//											InitialBurnedRadius,
											dipoleMoment, dipoleCenter, initBWithVectorPotential,
											&MetaData, triSphere);
	}

	TRACEF("HIERARCHY BUILT");

	if(projectChildrenToParents || 1)
	{
		LevelArrayIterator it = LevelArrayIterator(rhit.levelArray);
		bool usingDirectInit, usingVectorPotentialInit;

		{
			// A dummy call to get the values of the using direct/vector potential flags.
			grid *g = it.firstFromTop();
			g->MHDProfileInitializeGrid_B_and_CT_Fields(
					NULL, BurningTemperature, InitialBurnedRadius, dipoleMoment, dipoleCenter,
//						InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel) ? 1 : 0,
					initBWithVectorPotential, &MetaData, triSphere, &usingDirectInit, &usingVectorPotentialInit);
		}

		if(usingVectorPotentialInit)
		{
			for(grid* g = it.firstFromTop(); g; g = it.next())
				g->MHDProfileInitializeGrid_B_and_CT_Fields(
						&radialProfile, BurningTemperature, InitialBurnedRadius, dipoleMoment, dipoleCenter,
//						InitialBurnedRadius * (rhit.currentLevel == MaximumRefinementLevel) ? 1 : 0,
						initBWithVectorPotential, &MetaData, triSphere, NULL, NULL);
		}

		// Restore the consistency among levels by projecting each layer to its parent,
		// starting from the finest level.
		// We want to project the vector potential and take the curl afterwards.
		// Store the original project flags for the MHD fields and set them so
		// that the vector potential gets projected.
		it.projectChildrenToParents(usingDirectInit, usingVectorPotentialInit);

		if(usingVectorPotentialInit)
		{
			for(grid* g = it.firstFromTop(); g; g = it.next())
				g->MHDProfileInitializeGrid_CurlAndCenter(&radialProfile, BurningTemperature, InitialBurnedRadius,
															dipoleMoment, dipoleCenter, InitBWithVectorPotential,
															&MetaData, triSphere, usingDirectInit, usingVectorPotentialInit);
		}

		for(grid* g = it.firstFromTop(); g; g = it.next())
			g->MHDProfileInitializeGrid_TotalE_GasE(&radialProfile, BurningTemperature, InitialBurnedRadius, dipoleMoment,
											dipoleCenter, InitBWithVectorPotential, &MetaData, triSphere);

		for(grid* g = it.firstFromTop(); g; g = it.next())
			g->MHDProfileInitializeGrid5(&radialProfile, BurningTemperature, InitialBurnedRadius, dipoleMoment,
											dipoleCenter, InitBWithVectorPotential, &MetaData, triSphere);
	}

	if(triSphere)
	{
		// Free the tirSphere.
		triSphere->~TriSphere();
		delete triSphere;
		triSphere = NULL;
		// Turn off the static refinement shells
		StaticRefineShellInnerRadius[0] = -1;
		StaticRefineShellOuterRadius[0] = -2;
	}

	radialProfile.free();
	TRACEF("  Process initialization finished.");
	return SUCCESS;
}

int MHDProfileInitializeRestart(TopGridData * MetaData, LevelHierarchyEntry * *LevelArray, HierarchyEntry * TopGrid)
{
	// For MHD_RK and HD_RK
	NSpecies = 0;
	NColor = 1;
	NoMultiSpeciesButColors = 1;

	SphericalGravityComputePotential(LevelArray, MetaData, true);

	if(PerturbationOnRestart)
	{
		if(PerturbationMethod != 4)
		{
			ENZO_VFAIL("Perturbation-On-Restart is implemented for PerturbationMethod==4 only, but %lld found.",
						PerturbationMethod);
		}

		LevelArrayIterator it = LevelArrayIterator(LevelArray);
		TriSphere *triSphere = NewTriSphereFromParams(it.countLevels());
		for(grid* g = it.firstFromTop(); g; g = it.next())
		{
			g->PerturbWithTriSPhere(triSphere, NULL);
		}

		//Turn off perturbation for future restarts.
		PerturbationOnRestart = 0;
	}

	return SUCCESS;
}
