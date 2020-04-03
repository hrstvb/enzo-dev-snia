/*
 * MHDMagneticVectorPotentialInitialize.C
 *
 *  Created on: Jan 13, 2019
 *      Author: B. Hristov, D.C. Collins
 *
 *      A number of methods to set a vector potential
 *      for the initial magnetic field.  The vector potential is
 *      assigned to ElectricField.  One needs to call additionally
 *      this->MHD_Curl(...) and this->CenterMagneticField(...)
 *      to initialize the magnetic field.
 *
 *		When factor==0 the ElectricField is overwritten.  When
 *		factor!=0, the new vector potential is multiplied by
 *		*factor* and added to ElectricField.
 */

#include "math.h"
#include "myenzoutils.h"

#include "ErrorExceptions.h"
//#include "macros_and_parameters.h"
//#include "typedefs.h"
#include "global_data.h"
//#include "Fluxes.h"
//#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
//#include "Hierarchy.h"
//#include "TopGridData.h"
//#include "phys_constants.h"

int FindField(int f, int farray[], int n);

int grid::MHDClear_B_and_CT_Fields()
{
	if(!UseMHD)
		return SUCCESS;

	const int BFieldTypes[3] = { Bfield1, Bfield2, Bfield3 };
	int BFieldNum;
	float *field;

	for(int dim = 0; dim < 3; dim++)
	{
		if(0 < (BFieldNum = FindField(BFieldTypes[dim], FieldType, NumberOfBaryonFields)))
			if((field = BaryonField[BFieldNum]))
				arr_set(field, GetGridSize(), 0);

		if((field = MagneticField[dim]))
			arr_set(field, MagneticSize[dim], 0);

		if((field = ElectricField[dim]))
			arr_set(field, ElectricSize[dim], 0);
	}
	return SUCCESS;
}

int grid::InitializeMagneticUniformFieldVectorPotential(const float constMagneticField[3], const long float factor)
{
	const bool addToElectricField = true;
	const float f = ((addToElectricField) ? factor : 1);
	const float Bx = f * constMagneticField[0];
	const float By = f * constMagneticField[1];
	const float Bz = f * constMagneticField[2];

	bool has_j = GridRank > 1;
	bool has_k = GridRank > 2;
	size_t ni, nj, nk;
	float* EField = NULL;
	FLOAT Ex, Ey, Ez;

	int dim = 0;
	EField = ElectricField[dim];
	ni = ElectricDims[dim][0];
	nj = (has_j) ? ElectricDims[dim][1] : 1;
	nk = (has_k) ? ElectricDims[dim][2] : 1;
	if(has_k && By)
	{
		for(size_t k = 0; k < nk; k++)
		{
			Ex = - By * CellLeftEdge[2][k];
			for(size_t j = 0; j < nj; j++)
			{
				for(size_t i = 0; i < ni; i++)
				{
					if(addToElectricField)
						*EField++ += Ex;
					else
						*EField++ = Ex;
				}
			}
		}
	}
	else
	{
		if(!addToElectricField)
			arr_set(EField, ni * nj * nk, 0);
	}

	dim = 1;
	EField = ElectricField[dim];
	ni = ElectricDims[dim][0];
	nj = (has_j) ? ElectricDims[dim][1] : 1;
	nk = (has_k) ? ElectricDims[dim][2] : 1;
	if(has_j && Bz)
	{
		for(size_t k = 0; k < nk; k++)
		{
			for(size_t j = 0; j < nj; j++)
			{
				for(size_t i = 0; i < ni; i++)
				{
					Ey = Bz * CellLeftEdge[0][i];
					if(addToElectricField)
						*EField++ += Ey;
					else
						*EField++ = Ey;
				}
			}
		}
	}
	else
	{
		if(!addToElectricField)
			arr_set(EField, ni * nj * nk, 0);
	}

	dim = 2;
	EField = ElectricField[dim];
	ni = ElectricDims[dim][0];
	nj = (has_j) ? ElectricDims[dim][1] : 1;
	nk = (has_k) ? ElectricDims[dim][2] : 1;
	if(has_j && Bx)
	{
		for(size_t k = 0; k < nk; k++)
		{
			for(size_t j = 0; j < nj; j++)
			{
				Ez = Bx * CellLeftEdge[1][j];
				for(size_t i = 0; i < ni; i++)
				{
					if(addToElectricField)
						*EField++ += Ez;
					else
						*EField++ = Ez;
				}
			}
		}
	}
	else
	{
		if(!addToElectricField)
			arr_set(EField, ni * nj * nk, 0);
	}

	return SUCCESS;
}

int grid::InitializeMagneticUniformField(const float constMagneticField[3], const long float factor)
{
	const bool addToElectricField = true;
	const float f = ((addToElectricField) ? factor : 1);
	float b;

	float *BFields[3];
	MHD_SNIA_GetFields(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,BFields,NULL,NULL,NULL);

	for(int dim = 0; dim < 3; dim++)
	{
		if((b=f*constMagneticField[dim])==0)
			continue;

		arr_axpb(BFields[dim], gridSize, 1, b);
		arr_axpb(MagneticField[dim], MagneticSize[dim], 1, b);
	}

	return SUCCESS;
}

int grid::InitializeMagneticDipoleVectorPotential(const float dipoleMoment[3], const float dipoleCenter[3],
	const long float factor)
{
	const bool addToElectricField = true;
	const double factor_over_4pi = ((addToElectricField) ? factor : 1); // / (16 * atanl(1));
	const float mx = factor_over_4pi * dipoleMoment[0];
	const float my = factor_over_4pi * dipoleMoment[1];
	const float mz = factor_over_4pi * dipoleMoment[2];

	bool has_j = GridRank > 1;
	bool has_k = GridRank > 2;
	bool isCentered_i, isCentered_j, isCentered_k;
	size_t ni, nj, nk, gridSize, gridIndex;
	float* EField = NULL;
	FLOAT x, y, z;
	long double r, rx, ry, rz, rsquared, coeff, rz2, ry2_rz2, mx_ry, mx_rz, my_rz, my_rz__mz_ry, vp;
	long double minvp = 1e99, maxvp = 0;

	for(int dim = 0; dim < 3; dim++)
	{
		switch(dim){

		}

		ni = ElectricDims[dim][0];
		nj = (has_j) ? ElectricDims[dim][1] : 1;
		nk = (has_k) ? ElectricDims[dim][2] : 1;
		isCentered_i = ElectricDims[dim][0] == GridDimension[0];
		isCentered_j = ElectricDims[dim][1] == GridDimension[1];
		isCentered_k = ElectricDims[dim][2] == GridDimension[2];
		for(size_t k = 0; k < nk; k++)
		{
			z = (has_k) ? ((isCentered_k) ? CELLCENTER(2, k) : CellLeftEdge[2][k]) : 0;
			rz = z - dipoleCenter[2];
			rz2 = square(rz);
			switch(dim)
			{
			case 0:
				my_rz = my * rz;
				break;
			case 1:
				mx_rz = mx * rz;
				break;
			case 2:
				break;
			}

			for(size_t j = 0; j < nj; j++)
			{
				y = (has_j) ? ((isCentered_j) ? CELLCENTER(1, j) : CellLeftEdge[1][j]) : 0;
				ry = y - dipoleCenter[1];
				ry2_rz2 = square(ry) + rz2;
				switch(dim)
				{
				case 0:
					my_rz__mz_ry = my_rz - mz * ry;
					break;
				case 1:
					break;
				case 2:
					mx_ry = mx * ry;
					break;
				}

				gridIndex = ni * (j + nj * k);
				EField = ElectricField[dim] + gridIndex;
				for(size_t i = 0; i < ni; i++)
				{
					x = (isCentered_i) ? CELLCENTER(0, i) : CellLeftEdge[0][i];
					rx = x - dipoleCenter[0];
					rsquared = square(rx) + ry2_rz2;
					if(rsquared)
					{
						switch(dim)
						{
						case 0:
							vp = my_rz__mz_ry;
							break;
						case 1:
							vp = mz * rx - mx_rz;
							break;
						case 2:
							vp = mx_ry - my * rx;
							break;
						}
						vp /= rsquared;
					}
					else
					{
						vp = 0;
					}

					if(addToElectricField)
						*EField += vp;
					else
						*EField = vp;

					EField++;

				} //for i
			} //for j
		} //for k
	} //for dim

	return SUCCESS;
}

