/***********t************************************************************
 /
 /  GRID CLASS (COMPUTE TIME STEP)
 /
 /  written by: Greg Bryan
 /  date:       November, 1994
 /  modified1: 2010 Tom Abel, added MHD part
 /  modified2: 2015 Boyan Hristov added energy growth and burning fraction growth.
 /
 /  PURPOSE:
 /
 /  RETURNS:
 /    dt   - timestep
 /
 ************************************************************************/

// Compute the timestep from all the constrains for this grid.
//
// Somebody fix the error handling in this routine! please.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "myenzoutils.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferParameters.h"
#include "hydro_rk/EOS.h"
#include "hydro_rk/tools.h"
#include "phys_constants.h"
#include "LimitTimeStep.h"

/* function prototypes */

int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits, float *TimeUnits, float *VelocityUnits,
FLOAT Time);

extern "C" void PFORTRAN_NAME(calc_dt)(int *rank, int *idim, int *jdim, int *kdim, int *i1, int *i2, int *j1, int *j2,
int *k1, int *k2, hydro_method *ihydro, float *C2, FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
float *vgz, float *gamma, int *ipfree, float *aye, float *d, float *p, float *u, float *v, float *w, float *dt,
float *dtviscous, size_t *dtindex, size_t *dtviscousindex);

extern "C" void
FORTRAN_NAME(mhd_dt)(float *bxc, float *byc, float *bzc, float *vx, float *vy, float *vz, float *d, float *p,
float *gamma, float *dt, FLOAT *dx, FLOAT *dy, FLOAT *dz, int *idim, int *jdim, int *kdim, int * rank, int *i1,
int *i2, int *j1, int *j2, int *k1, int *k2, float* eng, size_t* dtindex);

float grid::populateDtLimitInfo(DtLimitInfo* info_inout, float dt)
{
	if(info_inout)
	{
		info_inout->dt = dt;
		info_inout->Grid = this;
		info_inout->GridID = this->ID;
		info_inout->ProcessorNumber = MyProcessorNumber;
		if(info_inout->dtIndex >= 0)
		{
			if(400 <= info_inout->reason)
			{
				getParticlePosition(info_inout->xyz, info_inout->dtIndex);
				//getParticleMass(info_inout->xyz, info_inout->dtIndex);
				//getParticleVelovcity(info_inout->xyz, info_inout->dtIndex);
				//etc.
			}
			else
			{
				get_ijk_xyz(info_inout->ijk, info_inout->xyz, info_inout->dtIndex);
				//get_rho_uvw
				//etc
			}
		}
	}

	return dt;
}

float grid::populateDtLimitInfo(DtLimitInfo* info_inout, float dt, dt_limit_reason reason, size_t dtIndex)
{
	if(info_inout)
	{
		info_inout->dt = dt;
		info_inout->dtIndex = dtIndex;
		info_inout->reason = reason;
	}

	return populateDtLimitInfo(info_inout, dt);
}

float grid::ComputeTimeStep()
{
	return ComputeTimeStep(NULL);
}

float grid::ComputeTimeStep(DtLimitInfo* dtLimitInfo_out)
{
	/* Return if this doesn't concern us. */

	if(ProcessorNumber != MyProcessorNumber)
		return populateDtLimitInfo(dtLimitInfo_out, huge_number, MAX_DT_NULL_REASON, 0);

	this->DebugCheck("ComputeTimeStep");

	/* initialize */

	float dt = huge_number;
	float dtTemp = dt;
	size_t dtIndex = MAX_DT_NO_INDEX_INFO;
	dt_limit_reason reason = MAX_DT_NULL_REASON;
	DtLimitInfo dtLimitInfo;
	dtLimitInfo.dt = dt;
	dtLimitInfo.dtIndex = dtIndex;
	dtLimitInfo.reason = reason;

	float dtAcceleration = huge_number;
	float dtHydro = huge_number;
	float dtHydroRK = huge_number;
	float dtBaryons = huge_number;
	float dtBurningEnergy = huge_number;
	float dtConduction = huge_number;
	float dtCooling = huge_number;
	float dtCR = huge_number;
	float dtExpansion = huge_number;
	float dtFreeFall = huge_number; //ProblemType==63 only
	float dtGasDrag = huge_number;
	float dtInternalEnergyGrowth = huge_number;
	float dtMHDLi = huge_number;
	float dtMHDRK = huge_number;
	float dtMHD = huge_number;
	float dtParticles = huge_number;
	float dtRadPressure = huge_number;
	float dtRadTransfer = huge_number;
	float dtSafetyVelocity = huge_number;
	float dtTotalEnergyGrowth = huge_number;
	float dtViscous = huge_number;
	size_t dtAccelerationIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtHydroIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtHydroRKIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtBaryonsIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtBurningEnergyIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtConductionIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtCoolingIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtCRIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtExpansionIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtFreeFallIndex = MAX_DT_NO_INDEX_INFO; //ProblemType==63 only
	size_t dtGasDragIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtInternalEnergyGrowthIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtMHDLiIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtMHDRKIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtMHDIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtParticlesIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtRadPressureIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtRadTransferIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtSafetyVelocityIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtTotalEnergyGrowthIndex = MAX_DT_NO_INDEX_INFO;
	size_t dtViscousIndex = MAX_DT_NO_INDEX_INFO;

	int dim, i, j, k, index, result;
	/* Compute the field size. */

	int size = GetGridSize();

	/* If using comoving coordinates, compute the expansion factor a.  Otherwise,
	 set it to one. */

	FLOAT a = 1, dadt;
	if(ComovingCoordinates)
		CosmologyComputeExpansionFactor(Time, &a, &dadt);
	float afloat = float(a);

	/* 1) Compute Courant condition for baryons. */

	if(NumberOfBaryonFields > 0 && (HydroMethod != HD_RK) && (HydroMethod != MHD_RK))
	{

		/* Find fields: density, total energy, velocity1-3. */

		int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum;
		this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum);

		/* For one-zone free-fall test, just compute free-fall time. */
		if(ProblemType == 63)
		{
			float *force_factor = new float[size];
			if(this->ComputeOneZoneCollapseFactor(force_factor) == FAIL)
			{
				ENZO_FAIL("Error in ComputeOneZoneCollapseFactor.\n");
			}

			dt = huge_number;
			for(k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
			{ // nothing
				for(j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
				{ // metallicity
					for(i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
					{ // energy

						index = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];

						dtTemp =
								POW(((3 * pi) / (32 * GravitationalConstant * BaryonField[DensNum][index] * (1 - force_factor[index]))),
									0.5);
						LimitDt(&dt, dtTemp, &dtFreeFallIndex, index);
					}
				}
			}

			delete[] force_factor;
			dtFreeFall = dt *= TestProblemData.OneZoneFreefallTimestepFraction;
			return populateDtLimitInfo(dtLimitInfo_out, dtFreeFall, MAX_DT_FREE_FALL_LIMITED, dtFreeFallIndex);
		}

		/* Compute the pressure. */

		float *pressure_field = new float[size];
		this->ComputePressure(Time, pressure_field, 0, 1); // Note: Force use of CRs to get sound speed correct

#ifdef UNUSED
		int Zero[3] =
		{	0,0,0}, TempInt[3] =
		{	0,0,0};
		for (dim = 0; dim < GridRank; dim++)
		TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */

		/* Call fortran routine to do calculation. */

		if(HydroMethod != MHD_Li)
		{
			PFORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension + 1, GridDimension + 2, GridStartIndex,
									GridEndIndex, GridStartIndex + 1, GridEndIndex + 1, GridStartIndex + 2,
									GridEndIndex + 2, &HydroMethod, &ZEUSQuadraticArtificialViscosity, CellWidth[0],
									CellWidth[1], CellWidth[2], GridVelocity, GridVelocity + 1, GridVelocity + 2,
									&Gamma, &PressureFree, &afloat, BaryonField[DensNum], pressure_field,
									BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], &dtHydro,
									&dtViscous, &dtHydroIndex, &dtViscousIndex);

			/* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
			dtHydro *= CourantSafetyNumber;
		}
		if(HydroMethod == MHD_Li)
		{
			/* 1.5) Calculate minimum dt due to MHD: Maximum Fast MagnetoSonic Shock Speed */

			//Cosmos nees this, for some reason.
			if(GridRank < 3)
			{
				if(CellWidth[2] == NULL)
					CellWidth[2] = new FLOAT;
				CellWidth[2][0] = 1.0;
				if(GridRank < 2)
				{
					if(CellWidth[1] == NULL)
						CellWidth[1] = new FLOAT;
					CellWidth[1][0] = 1.0;
				}
			}
			int Rank_Hack = 3; //MHD needs a 3d timestep always.
			FORTRAN_NAME(mhd_dt)(BaryonField[B1Num], BaryonField[B2Num], BaryonField[B3Num], BaryonField[Vel1Num],
									BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[DensNum], pressure_field,
									&Gamma, &dtMHDLi, CellWidth[0], CellWidth[1], CellWidth[2], GridDimension,
									GridDimension + 1, GridDimension + 2, &Rank_Hack, GridStartIndex, GridEndIndex,
									GridStartIndex + 1, GridEndIndex + 1, GridStartIndex + 2, GridEndIndex + 2,
									BaryonField[TENum], &dtMHDLiIndex);

			dtMHDLi *= CourantSafetyNumber;
			dtMHDLi *= afloat;
		} //if HydroMethod== MHD_Li

		/* Clean up */
		delete[] pressure_field;
	}

	if(NumberOfBaryonFields > 0 && (HydroMethod == HD_RK))
	{

		int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
		if(this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum) == FAIL)
		{
			fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
			exit(FAIL);
		}

		FLOAT dxinv = 1.0 / CellWidth[0][0] / a;
		FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] / a : 0.0;
		FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] / a : 0.0;
		float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
		float rho, p, vx, vy, vz, v2, eint, etot, h, cs, dpdrho, dpde, v_signal_x, v_signal_y, v_signal_z;
		int n = 0;
		for(k = 0; k < GridDimension[2]; k++)
		{
			for(j = 0; j < GridDimension[1]; j++)
			{
				for(i = 0; i < GridDimension[0]; i++, n++)
				{
					rho = BaryonField[DensNum][n];
					vx = BaryonField[Vel1Num][n];
					vy = BaryonField[Vel2Num][n];
					vz = BaryonField[Vel3Num][n];

					if(DualEnergyFormalism)
					{
						eint = BaryonField[GENum][n];
					}
					else
					{
						etot = BaryonField[TENum][n];
						v2 = vx * vx + vy * vy + vz * vz;
						eint = etot - 0.5 * v2;
					}

					EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);

					v_signal_y = v_signal_z = 0;

					v_signal_x = (cs + fabs(vx));
					if(GridRank > 1)
						v_signal_y = (cs + fabs(vy));
					if(GridRank > 2)
						v_signal_z = (cs + fabs(vz));

					dt_x = v_signal_x * dxinv;
					dt_y = v_signal_y * dyinv;
					dt_z = v_signal_z * dzinv;

					dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

					if(dt_ltemp > dt_temp)
					{
						dt_temp = dt_ltemp;
						dtHydroRKIndex = ELT(i, j, k);
					}
				}
			}
		}

		dtHydroRK = CourantSafetyNumber / dt_temp;

	}

	// MHD
	if(NumberOfBaryonFields > 0 && HydroMethod == MHD_RK)
	{

		int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum;
		if(this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num,
											PhiNum) == FAIL)
			ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

		FLOAT dxinv = 1.0 / CellWidth[0][0] / a;
		FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] / a : 0.0;
		FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] / a : 0.0;
		float vxm, vym, vzm, Bm, rhom;
		float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
		float rho, p, vx, vy, vz, v2, eint, etot, h, cs, cs2, dpdrho, dpde, v_signal_x, v_signal_y, v_signal_z, cf, cf2,
				temp1, Bx, By, Bz, B2, ca2;
		int n = 0;
		float rho_dt, B_dt, v_dt;
		for(k = 0; k < GridDimension[2]; k++)
		{
			for(j = 0; j < GridDimension[1]; j++)
			{
				for(i = 0; i < GridDimension[0]; i++, n++)
				{
					rho = BaryonField[DensNum][n];
					vx = BaryonField[Vel1Num][n];
					vy = BaryonField[Vel2Num][n];
					vz = BaryonField[Vel3Num][n];
					Bx = BaryonField[B1Num][n];
					By = BaryonField[B2Num][n];
					Bz = BaryonField[B3Num][n];

					B2 = Bx * Bx + By * By + Bz * Bz;
					if(DualEnergyFormalism)
					{
						eint = BaryonField[GENum][n];
					}
					else
					{
						etot = BaryonField[TENum][n];
						v2 = vx * vx + vy * vy + vz * vz;
						eint = etot - 0.5 * v2 - 0.5 * B2 / rho;
					}

					v_signal_y = v_signal_z = 0;

					EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
					cs2 = cs * cs;
					temp1 = cs2 + B2 / rho;

					ca2 = Bx * Bx / rho;
					cf2 = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4.0 * cs2 * ca2));
					cf = sqrt(cf2);
					v_signal_x = (cf + fabs(vx));

					if(GridRank > 1)
					{
						ca2 = By * By / rho;
						cf2 = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4.0 * cs2 * ca2));
						cf = sqrt(cf2);
						v_signal_y = (cf + fabs(vy));
					}

					if(GridRank > 2)
					{
						ca2 = Bz * Bz / rho;
						cf2 = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4.0 * cs2 * ca2));
						cf = sqrt(cf2);
						v_signal_z = (cf + fabs(vz));
					}

					dt_x = v_signal_x * dxinv;
					dt_y = v_signal_y * dyinv;
					dt_z = v_signal_z * dzinv;

					dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

					if(dt_ltemp > dt_temp)
					{
						dt_temp = dt_ltemp;
						rho_dt = rho;
						B_dt = sqrt(Bx * Bx + By * By + Bz * Bz);
						v_dt = max(fabs(vx), fabs(vy));
						v_dt = max(v_dt, fabs(vz));
						dtMHDRKIndex = ELT(i, j, k);
					}

				}
			}
		}
		dtMHDRK = CourantSafetyNumber / dt_temp;
		//    fprintf(stderr, "ok %g %g %g\n", dt_x,dt_y,dt_z);
		//    if (dtMHD*TimeUnits/yr < 5) {
		//float ca = B_dt/sqrt(rho_dt)*VelocityUnits;
		//printf("dt=%g, rho=%g, B=%g\n, v=%g, ca=%g, dt=%g", dtMHD*TimeUnits/yr, rho_dt*DensityUnits, B_dt*MagneticUnits,
		//    v_dt*VelocityUnits/1e5, ca/1e5, LengthUnits/(dxinv*ca)/yr*CourantSafetyNumber);
		// }

	} // HydroMethod = MHD_RK

	/* 2) Calculate dt from particles. */

	if(NumberOfParticles > 0)
	{

		/* Compute dt constraint from particle velocities. */

		for(dim = 0; dim < GridRank; dim++)
		{
			float dCell = CellWidth[dim][0] * a;
			for(i = 0; i < NumberOfParticles; i++)
			{
				dtTemp = dCell / max(fabs(ParticleVelocity[dim][i]), tiny_number);
//	dtParticles = min(dtParticles, dtTemp);
				if(dtParticles > dtTemp)
				{
					dtParticles = dtTemp;
					dtParticlesIndex = i;
				}
			}
		}

		/* Multiply resulting dt by ParticleCourantSafetyNumber. */

		dtParticles *= ParticleCourantSafetyNumber;

	}

	/* 3) Find dt from expansion. */

	if(ComovingCoordinates)
	{
		if(CosmologyComputeExpansionTimestep(Time, &dtExpansion) == FAIL)
		{
			fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
			exit(FAIL);
		}
		dtExpansionIndex = MAX_DT_INDEX_NOT_APPLICABLE; //Not applicable
	}
	/* 4) Calculate minimum dt due to acceleration field (if present). */

	if(SelfGravity)
	{
		for(dim = 0; dim < GridRank; dim++)
			if(AccelerationField[dim])
				for(i = 0; i < size; i++)
				{
					dtTemp = sqrt(CellWidth[dim][0] / fabs(AccelerationField[dim][i]) + tiny_number);
					//	      dtAcceleration = min(dtAcceleration, dtTemp);
					LimitDt(&dtAcceleration, dtTemp, &dtAccelerationIndex, i);
				}

		if(dtAcceleration != huge_number)
			dtAcceleration *= 0.5;
	}

	/* 5) Calculate minimum dt due to thermal conduction. */

	if(IsotropicConduction || AnisotropicConduction)
	{
		if(this->ComputeConductionTimeStep(dtConduction) == FAIL)
			ENZO_FAIL("Error in ComputeConductionTimeStep.\n");

		dtConduction *= float(NumberOfGhostZones);     // for subcycling
		dtConductionIndex = MAX_DT_NO_INDEX_INFO; //TODO
	}

	/* 6) Calculate minimum dt due to CR diffusion */

	if(CRModel && CRDiffusion)
	{
		if(this->ComputeCRDiffusionTimeStep(dtCR) == FAIL)
		{
			fprintf(stderr, "Error in ComputeCRDiffusionTimeStep.\n");
			return FAIL;
		}
		dtCR *= CRCourantSafetyNumber;
		dtCR *= float(NumberOfGhostZones);  // for subcycling
		dtCRIndex = MAX_DT_NO_INDEX_INFO; //TODO
	}

	/* 7) GasDrag time step */
	if(UseGasDrag && GasDragCoefficient != 0.)
	{
		dtGasDrag = 0.5 / GasDragCoefficient;
		dtGasDragIndex = MAX_DT_INDEX_NOT_APPLICABLE; //Not applicable
	}

	/* Cooling time */
	if(UseCoolingTimestep == TRUE)
	{
		float *cooling_time = new float[size];
		if(this->ComputeCoolingTime(cooling_time, TRUE) == FAIL)
		{
			ENZO_FAIL("Error in grid->ComputeCoolingTime.\n");
		}

		for(k = GridStartIndex[2]; k < GridEndIndex[2]; k++)
		{
			for(j = GridStartIndex[1]; j < GridEndIndex[1]; j++)
			{
				index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
				for(i = GridStartIndex[0]; i < GridEndIndex[0]; i++, index++)
					LimitDt(&dtCooling, cooling_time[index], &dtCoolingIndex, index);
			}
		}
		dtCooling *= CoolingTimestepSafetyFactor;

		delete[] cooling_time;
	}

	/* 8) calculate minimum timestep */

//	dt = min(    dtBaryons, dtParticles   );
//	dt = min(dt, dtMHD         );
//	dt = min(dt, dtViscous     );
//	dt = min(dt, dtAcceleration);
//	dt = min(dt, dtExpansion   );
//	dt = min(dt, dtConduction  );
//	dt = min(dt, dtCR          );
//	dt = min(dt, dtGasDrag     );
//	dt = min(dt, dtCooling     );
#ifdef TRANSFER

	float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, TimeUnits, aUnits = 1;

	if(GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time) == FAIL)
	{
		ENZO_FAIL("Error in GetUnits.");
	}

	/* 8) If using radiation pressure, calculate minimum dt */

	float absVel, absAccel;

	if(RadiationPressure && RadiativeTransfer)
	{

		int RPresNum1, RPresNum2, RPresNum3;
		if(IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) == FAIL)
		{
			ENZO_FAIL("Error in IdentifyRadiationPressureFields.");
		}

		for(i = 0; i < size; i++)
			for(dim = 0; dim < GridRank; dim++)
			{
				dtTemp = sqrt(CellWidth[dim][0] / (fabs(BaryonField[RPresNum1 + dim][i]) + tiny_number));
//				dtRadPressure = min(dtRadPressure, dtTemp);
				LimitDt(&dtRadPressure, dtTemp, &dtRadPressureIndex, i);
			}

		if(dtRadPressure < huge_number)
			dtRadPressure *= 0.5;

//		dt = min(dt, dtRadPressure);

	} /* ENDIF RadiationPressure */

	/* 9) Safety Velocity to limit timesteps */
	if(TimestepSafetyVelocity > 0)
	{
		dtSafetyVelocity = a * CellWidth[0][0] / (TimestepSafetyVelocity * 1e5 / VelocityUnits);    // parameter in km/s
		dtSafetyVelocityIndex = MAX_DT_INDEX_NOT_APPLICABLE;
//	dt = min(dt, dtSafetyVelocity);
	}
	/* 10) FLD Radiative Transfer timestep limitation */
	if(RadiativeTransferFLD)
//		dt = min(dt, MaxRadiationDt);
		dtRadTransfer = MaxRadiationDt;
	dtRadTransferIndex = MAX_DT_INDEX_NOT_APPLICABLE;

#endif /* TRANSFER */

	/* 12) Compute time step from energy growth limits. */			//[BH]
//if( ComputeEnergyGrowthTimeStep( dEdtParent ) == FAIL )			//[BH]
//	ENZO_FAIL("Error in ComputeEnergyGrowthTimeStep.\n");			//[BH]
// dtEnergyGrowth is calculated somwhere else since				//[BH]
// OldBaryonField is not available at this point.				//[BH]
	/* 13) Compute time step from diffusion of burned fraction. *///[BH]
	if(UseBurning)
	{								//[BH]
		if(ComputeBurningFractionDiffusionTimeStep(&dtBurningEnergy) == FAIL)	//[BH]
			ENZO_FAIL("Error in ComputeBurnedFractionTimeStep.\n");		//[BH]
//  	dt = min(dt, dtBurnedFractioin);					//[BH]
		dtBurningEnergyIndex = MAX_DT_INDEX_NOT_APPLICABLE;
	}										//[BH]
	/* Debugging info. */

	LimitDt(&dt, dtAcceleration, &dtLimitInfo, MAX_DT_ACCELERATION_LIMITED, dtAccelerationIndex);
	LimitDt(&dt, dtHydro, &dtLimitInfo, MAX_DT_HYDRO_LIMITED, dtHydroIndex);
	LimitDt(&dt, dtHydroRK, &dtLimitInfo, MAX_DT_HYDRO_RK_LIMITED, dtHydroRKIndex);
	LimitDt(&dt, dtConduction, &dtLimitInfo, MAX_DT_CONDUCTION_LIMITED, dtConductionIndex);
	LimitDt(&dt, dtCooling, &dtLimitInfo, MAX_DT_COOLING_LIMITED, dtCoolingIndex);
	LimitDt(&dt, dtCR, &dtLimitInfo, MAX_DT_CR_LIMITED, dtCRIndex);
	LimitDt(&dt, dtExpansion, &dtLimitInfo, MAX_DT_EXPANSION_LIMITED, dtExpansionIndex);
	LimitDt(&dt, dtGasDrag, &dtLimitInfo, MAX_DT_GASDRAG_LIMITED, dtGasDragIndex);
	LimitDt(&dt, dtMHDLi, &dtLimitInfo, MAX_DT_MHD_LI_LIMITED, dtMHDLiIndex);
	LimitDt(&dt, dtMHDRK, &dtLimitInfo, MAX_DT_MHD_RK_LIMITED, dtMHDRKIndex);
	LimitDt(&dt, dtParticles, &dtLimitInfo, MAX_DT_PARTICLES_LIMITED, dtParticlesIndex);
	LimitDt(&dt, dtViscous, &dtLimitInfo, MAX_DT_VISCOUS_LIMITED, dtViscousIndex);
	LimitDt(&dt, dtRadPressure, &dtLimitInfo, MAX_DT_RAD_PRESSURE_LIMITED, dtRadPressureIndex);
	LimitDt(&dt, dtSafetyVelocity, &dtLimitInfo, MAX_DT_SAFETY_VELOCITY_LIMITED, dtSafetyVelocityIndex);
	LimitDt(&dt, dtRadTransfer, &dtLimitInfo, MAX_DT_RAD_TRANSFER_LIMITED, dtRadTransferIndex);
	LimitDt(&dt, dtTotalEnergyGrowth, &dtLimitInfo, MAX_DT_TOTAL_ENERGY_GROWTH_LIMITED, dtTotalEnergyGrowthIndex); //[BH]
	LimitDt(&dt, dtInternalEnergyGrowth, &dtLimitInfo, MAX_DT_INTERNAL_ENERGY_GROWTH_LIMITED,
			dtInternalEnergyGrowthIndex); //[BH]
	LimitDt(&dt, dtBurningEnergy, &dtLimitInfo, MAX_DT_BURNING_ENERGY_LIMITED, dtBurningEnergyIndex);

	if(debug1)
	{
		printf("ComputeTimeStep = %"ESYM" (", dt);
		if(HydroMethod != MHD_RK && HydroMethod != MHD_Li && NumberOfBaryonFields > 0)
		{
			printf("Bar = %"ESYM" ", dtHydro);
			printf("Hyd_RK = %"ESYM" ", dtHydroRK);
		}
		if(HydroMethod == MHD_Li)
			printf("dtMHD_Li = %"ESYM" ", dtMHDLi);
		if(HydroMethod == MHD_RK)
			printf("dtMHD_RK = %"ESYM" ", dtMHDRK);
		if(HydroMethod == Zeus_Hydro)
			printf("Vis = %"ESYM" ", dtViscous);
		if(ComovingCoordinates)
			printf("Exp = %"ESYM" ", dtExpansion);
		if(dtAcceleration != huge_number)
			printf("Acc = %"ESYM" ", dtAcceleration);
		if(NumberOfParticles)
			printf("Part = %"ESYM" ", dtParticles);
		if(UseCoolingTimestep)
			printf("Cool = %"ESYM" ", dtCooling);
		if(IsotropicConduction || AnisotropicConduction)
			printf("Cond = %"ESYM" ", (dtConduction));
		if(UseGasDrag)
			printf("Drag = %"ESYM" ", (dtGasDrag));
		if(dtInternalEnergyGrowth < huge_number)			//[BH]
			printf("dE/E = %"ESYM" ", dtInternalEnergyGrowth);	//[BH]
		if(dtTotalEnergyGrowth < huge_number)			//[BH]
			printf("dE/E = %"ESYM" ", dtTotalEnergyGrowth);	//[BH]
		//TODO burning					//[BH]
		//if( UseBurning )					//[BH]
		//  printf("df = %"ESYM" ", dtBurnedFraction);	//[BH]
		printf(")\n");
	}

	*dtLimitInfo_out = dtLimitInfo;
	return populateDtLimitInfo(dtLimitInfo_out, dt);
}

int grid::ComputeEnergyGrowthTimeStep()
{
	ComputeTotalEnergyGrowthTimeStep(&dtTotalEnergyGrowth, &dtTotalEnergyGrowthIndex);
	ComputeInternalEnergyGrowthTimeStep(&dtInternalEnergyGrowth, &dtInternalEnergyGrowthIndex);
	return SUCCESS;
}

int grid::ComputeTotalEnergyGrowthTimeStep(float* dt, size_t* dtIndex)
{
	if(TotalEnergyRelativeGrowthLimit <= 0)
		return SUCCESS;

	if(dtFixed <= tiny_number)
		return SUCCESS;

	int TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
	if(TENum < 0)
		ENZO_FAIL("Grid::ComputeEnergyGrowthTimeStep: Cannot find total energy.");

	float *oldTE = OldBaryonField[TENum];
	if(oldTE == NULL)
		return SUCCESS;

	int useGE = DualEnergyFormalism;
	float *GE, *oldGE;

	float *TE = BaryonField[TENum];
	float dE, dEMax = -1;

	size_t i1 = GridStartIndex[0];
	size_t i2 = GridEndIndex[0];
	size_t j1 = (GridRank > 1) ? 0 : GridStartIndex[1];
	size_t j2 = (GridRank > 1) ? 0 : GridEndIndex[1];
	size_t k1 = (GridRank > 2) ? 0 : GridStartIndex[2];
	size_t k2 = (GridRank > 2) ? 0 : GridEndIndex[2];
	size_t index = 0, maxi = 0;
	for(size_t k = k1; k <= k2; k++)
	{
		for(size_t j = j1; j <= j2; j++)
		{
			index = i1 + GridDimension[0] * (j + GridDimension[1] * k);
			for(size_t i = i1; i <= i2; i++)
			{
				if((dE = oldTE[i]) < tiny_number)
				{
					if(dEMax < (dE = abs(TE[i] / dE)))
					{
						dEMax = dE;
						maxi = i;
					}
				}
				index++;
			}
		}
	}

	*dt = TotalEnergyRelativeGrowthLimit * dtFixed / dEMax;
	*dtIndex = maxi;

	return SUCCESS;
} //end grid::ComputeEnergyGrowthTimeStep

int grid::ComputeInternalEnergyGrowthTimeStep(float* dt, size_t* dtIndex)
{
	if(TotalEnergyRelativeGrowthLimit <= 0)
		return SUCCESS;

	if(dtFixed <= tiny_number)
		return SUCCESS;

	int DensNum, GENum, TENum, VxNum, VyNum, VzNum, BxNum, ByNum, BzNum;
	float *GE, *TE, *Dens, *Vx, *Vy, *Vz, *Bx, *By, *Bz;
	float *oldGE, *oldTE, *oldDens, *oldVx, *oldVy, *oldVz, *oldBx, *oldBy, *oldBz;
	float gasE, oldGasE, dE, dEMax = 0;
	int hasGE = DualEnergyFormalism;
	bool hasMHD = UseMHD || UseMHDCT;

	if(hasGE)
	{
		GENum = FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
		if(!(GE = BaryonField[GENum]))
			ENZO_FAIL("grid::ComputeInternalEnergyGrowthTimeStep: Error: Can't find internal energy field.");
		if(!(oldGE = OldBaryonField[GENum]))
			return SUCCESS;
	}
	else
	{
		if(hasMHD)
		{
			if(IdentifyPhysicalQuantities(DensNum, GENum, VxNum, VyNum, VzNum, TENum, BxNum, ByNum, BzNum) == FAIL)
				ENZO_FAIL("Error in IdentifyPhysicalQuantities (with B-fields).");

			Bx = BaryonField[BxNum];
			By = BaryonField[ByNum];
			Bz = BaryonField[BzNum];
		}
		else
		{
			if(IdentifyPhysicalQuantities(DensNum, GENum, VxNum, VyNum, VzNum, TENum) == FAIL)
				ENZO_FAIL("Error in IdentifyPhysicalQuantities (no B-fields).");
		}

		Dens = BaryonField[DensNum];
		TE = BaryonField[TENum];
		Vx = BaryonField[VxNum];
		if(GridRank > 1)
			Vy = BaryonField[VyNum];
		if(GridRank > 2)
			Vz = BaryonField[VzNum];
	}

	size_t i1 = GridStartIndex[0];
	size_t i2 = GridEndIndex[0];
	size_t j1 = (GridRank > 1) ? 0 : GridStartIndex[1];
	size_t j2 = (GridRank > 1) ? 0 : GridEndIndex[1];
	size_t k1 = (GridRank > 2) ? 0 : GridStartIndex[2];
	size_t k2 = (GridRank > 2) ? 0 : GridEndIndex[2];
	size_t index = 0, maxi = 0;
	for(size_t k = k1; k <= k2; k++)
	{
		for(size_t j = j1; j <= j2; j++)
		{
			index = i1 + GridDimension[0] * (j + GridDimension[1] * k);
			for(size_t i = i1; i <= i2; i++, index++)
			{
				if(hasGE)
				{
					gasE = GE[index];
					oldGasE = oldGE[index];
				}
				else
				{
					gasE = TE[index];
					oldGasE = oldTE[index];
					if(hasMHD)
					{
						gasE -= (square(Bx[index]) + square(By[index]) + square(Bz[index])) / Dens[index];
						oldGasE -= (square(oldBx[index]) + square(oldBy[index]) + square(oldBz[index]))
								/ oldDens[index];
					}
					switch(GridRank)
					{
					case 3:
						gasE -= square(Vz[index]);
						oldGasE -= square(oldVz[index]);
						// no break, fall through
					case 2:
						gasE -= square(Vy[index]);
						oldGasE -= square(oldVy[index]);
						// no break, fall through
					}
					gasE -= square(Vx[index]);
					oldGasE -= square(oldVx[index]);
				}
				if(oldGasE > tiny_number)
				{
					if(dEMax < (dE = abs(gasE / oldGasE)))
					{
						dEMax = dE;
						maxi = i;
					}
				}
			}
		}
	}
	return SUCCESS;
}
