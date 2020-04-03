//
// ExtraFunction
// Generic grid member function for debugging.
//

#include <stdio.h>
#include <cstdarg>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "DebugMacros.h"

int grid::ExtraFunction(char * message)
{
	return SUCCESS;
}

long long grid::ExtraFunction(char *message, long long intarg, ...)
{
	/*
	 *	std::va_list args;
	 *	va_start(args, myCount);
	 *	for (int i = 0; i < myCount; ++i)
	 *	{
	 *		// Print integers
	 *		std::cout << va_arg(args, int) << '\n';
	 *	}
	 *	va_end(args);
	 */
	static int numCalls = 0;
	numCalls++;
	if(numCalls >= 21)
		return 0;

	size_t NI = GridDimension[0];
	size_t NJ = GridDimension[1];
	size_t NK = GridDimension[2];
	size_t i_first = 0;
	size_t i_last = NI - 1;
	size_t i_count = i_last - i_first + 1;
	size_t J = NJ / 2;
	size_t K = NK / 2;
	size_t index_first = ELT(i_first, J, K);
	FLOAT DX = CellWidth[0][0];
	size_t index;
	float *field = NULL;
	char *fieldName = NULL;
	int cycleNum = -1;
	FLOAT t = 0;
	int vaCount = 0;

	cycleNum = intarg;
//	std::va_list args;
//	va_start(args, vaCount);
//	t = va_arg(args, FLOAT);
//	va_end(args);

	float *rhoField, *totEField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;
	float *presField = NULL, *mgp = NULL, *mgpor = NULL, *g_mgpor = NULL;
	float *gxField = AccelerationField[0];
	float *gyField = AccelerationField[1];
	float *gzField = AccelerationField[2];

	MHD_SNIA_GetFields(&rhoField, &totEField, &gasEField, &vxField, &vyField, &vzField, NULL, &BxField, &ByField,
						&BzField, NULL, &rhoNiField, NULL, NULL);

	{
		FILE *outfile = stderr;
		char fname[128];
		sprintf(fname, "extrafunc%02lld.txt", 0);
		outfile = fopen(fname, "w");
		fprintf(outfile, "i j k x y z r");
		for(int i = 1; i <= numCalls; i++)
			fprintf(outfile, " g%02ldd -dP%02ldd (g+dP)%02ldd", i, i, i);
		fprintf(outfile, "\n");
		fclose(outfile);

		sprintf(fname, "extrafunc%02lld.txt", numCalls);
		TRACEGF(" %s %s cycle %lld callnum %lld", fname, message, cycleNum, numCalls);
		outfile = fopen(fname, "w");
		fprintf(outfile, "%s cycle %lld callnum %lld\n", message, cycleNum, numCalls);

		float rr[200], dd[200], gg[200], PP[200], dPdP[200];
		for(int i = 0; i < NI; i++)
		{
			int j = NJ / 2;
			int k = i;
			index = ELT(i, j, k);

			FLOAT x = CELLCENTER(0, i);
			FLOAT y = CELLCENTER(1, j);
			FLOAT z = CELLCENTER(2, k);
			FLOAT r = lenl(x, y, z);
			float g = 0, v = 0, b = 0, P = 0;
			if(gxField)
				g = lenl(gxField[index], gyField[index], gzField[index]);
			if(vxField==NULL)
				break;
			v = lenl(vxField[index], vyField[index], vzField[index]);
			if(BxField)
				b = lenl(BxField[index], ByField[index], BzField[index]);
			P = totEField[index] * rhoField[index] * (Gamma - 1);

			rr[i] = r;
			dd[i] = rhoField[index];
			gg[i] = g;
			PP[i] = P;
			dPdP[i] = 0;
		}

		for(int i = 1; i < NI - 1; i++)
		{
			int j = NJ / 2;
			int k = i;
			index = ELT(i, j, k);

			dPdP[i] = (PP[i + 1] - PP[i - 1]) / (rr[i + 1] - rr[i - 1]) / dd[i];
		}

		char tabs[90];
		arr_set(tabs, 90, '\t');
		tabs[3 * numCalls] = '\0';
		for(int i = 1; i < NI - 1; i++)
		{
			int j = NJ / 2;
			int k = i;
			index = ELT(i, j, k);

			FLOAT x = CELLCENTER(0, i);
			FLOAT y = CELLCENTER(1, j);
			FLOAT z = CELLCENTER(2, k);
			fprintf(outfile, "%lld %lld %lld   %e %e %e %e%s%e %e %e\n", i, j, k, x, y, z, rr[i], tabs, gg[i], -dPdP[i],
					gg[i] + dPdP[i]);
		}
//		fprintf(outfile, "i j k x y z r rho e P rhoni gx gy gz g vx vy vz v Bx By Bz B\n");
		//		for(int i = 0; i < NI; i++)
		//		{
		//			int j = NJ/2;
		//			int k = i;
		//			index = ELT(i, j, k);
		//
		//			FLOAT x = CELLCENTER(0, i);
		//			FLOAT y = CELLCENTER(1, j);
		//			FLOAT z = CELLCENTER(2, k);
		//			FLOAT r = lenl(x, y, z);
		//			float g = lenl(gxField[index], gyField[index], gzField[index]);
		//			float v = lenl(vxField[index], vyField[index], vzField[index]);
		//			float b = lenl(BxField[index], ByField[index], BzField[index]);
		//			float P = totEField[index] * rhoField[index] * (Gamma - 1);
		//			float dPdr
		//			fprintf(outfile, "%lld %lld %lld   %e %e %e %e    %e %e %e %e   %e %e %e %e   %e %e %e %e   %e %e %e %e\n",
		//				i, j, k,
		//				x,y,z,r,
		//				rhoField  [index],
		//				totEField [index],
		//				P,
		//				rhoNiField[index],
		//				gxField[index],
		//				gyField[index],
		//				gzField[index],
		//				g,
		//				vxField   [index],
		//				vyField   [index],
		//				vzField   [index],
		//				v,
		//				BxField   [index],
		//				ByField   [index],
		//				BzField   [index],
		//				b
		//				);
		//		}

		fclose(outfile);
		return 0;
	}

	switch(intarg)
	{
	case 0:
		field = vxField;
		fieldName = "V_x";
		break;
	case 1:
		field = AccelerationField[0];
		fieldName = "g_x";
		break;
	case 2:
		field = totEField;
		fieldName = "totE";
		break;
	case 3:
		fieldName = "-gradP/rho";
		break;
	case 4:
		fieldName = "(g-gradP/rho)";
		break;
	case 5:
		fieldName = "TotalBindingEnergy";
		break;
	}

	if(intarg == 5)
	{
		double totalBE = 0;
		for(int k = 0; k < GridDimension[2]; k++)
		{
			for(int j = 0; j < GridDimension[1]; j++)
			{
				for(int i = 0; i < GridDimension[0]; i++)
				{
					FLOAT r = lenl(CELLCENTER(0, i), CELLCENTER(1, j), CELLCENTER(2, k));
					if(r > OuterVelocitiesSphereRadius)
						continue;

					size_t index = ELT(i, j, k);
					double be = lenl(AccelerationField[0][index], AccelerationField[1][index],
										AccelerationField[2][index]);

					be *= r * rhoField[index] * getCellVolume();
					totalBE += be;
				}
			}
		}
		TRACEGF(" Total binding energy within radius  %e  at t=  %e  is  %e", OuterVelocitiesSphereRadius, t, totalBE);
		return 0;
	}

	if(intarg != 4)
		return 0;

	{
		field = new float[i_count];
		presField = new float[i_count];
		mgp = new float[i_count];
		mgpor = new float[i_count];
		g_mgpor = new float[i_count];

		index = index_first;
		for(size_t i = i_first; i <= i_last; i++)
		{
			presField[i] = rhoField[index] * totEField[index] * (Gamma - 1);
			index++;
		}

		index_first++;
		i_first++;
		i_last--;

		index = index_first;
		for(size_t i = i_first; i <= i_last; i++)
		{
			float P1 = presField[i - 1];
			float P2 = presField[i + 1];
			mgp[i] = -(P2 - P1) / 2 / DX; //Minus grap P
			mgpor[i] = mgp[i] / rhoField[index]; //Minus grad P over rho
			g_mgpor[i] = mgpor[i] + AccelerationField[0][index];
			index++;
		}
	}

	FILE *outFile = stderr;
	char *FMODE_W = "w+";
	char *FMODE_A = "a";
	char *fopenmode = (numCalls) ? FMODE_A : FMODE_W;
	outFile = fopen("extrafunc.txt", fopenmode);

	{
		const size_t MAXTABS = 256;
		char TABS[MAXTABS];
		for(int i = 0; i < MAXTABS; i++)
			TABS[i] = '\0';
		for(int i = 0; i < numCalls; i++)
			strcat(TABS, "\t\t\t\t\t\t");
		fprintf(outFile, "\t\tx"
				"\tg-dP/dx/r(cyc=0lev=0)\tg\tdP/dx/rho\tdP/dx\tP\trho"
				"\tg-dP/dx/r(cyc=0lev=1)\tg\tdP/dx/rho\tdP/dx\tP\trho"
				"\tg-dP/dx/r(cyc=1lev=0)\tg\tdP/dx/rho\tdP/dx\tP\trho"
				"\tg-dP/dx/r(cyc=1lev=1)\tg\tdP/dx/rho\tdP/dx\tP\trho\n");
		fprintf(outFile, "\n%s__callNum_%lld\n", message, numCalls);
		fprintf(outFile, "%s__i__x__a__g__dpdxrho__dpdx__p__rho=\n", message);
		index = index_first;
		for(int i = i_first; i <= i_last; i++)
		{
			fprintf(outFile, "%s:\t%lld\t%e\t%s%e\t%e\t%e\t%e\t%e\t%e\n", message, i, CELLCENTER(0, i), TABS,
					g_mgpor[i], AccelerationField[0][index], mgpor[i], mgp[i], presField[i], rhoField[index]);
			index++;
		}
		fprintf(outFile, "\n");
	}

	fclose(outFile);

	return 0;
}
