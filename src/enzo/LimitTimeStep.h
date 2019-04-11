/*
 * This file includes utilities used in computing the time step
 * for the grids, levels and hierarchy. It helps track the reason
 * and the location of the time step limit.
 */
#ifndef LIMIT_TIME_STEP_H
#define LIMIT_TIME_STEP_H

#include "stddef.h"
#include "macros_and_parameters.h"
#include "global_data.h"

struct HierarchyEntry;
class grid;

const size_t MAX_DT_NO_INDEX_INFO = -2;
const size_t MAX_DT_INDEX_NOT_APPLICABLE = -1;

typedef int dt_limit_reason;
const dt_limit_reason MAX_DT_NULL_REASON = 0;
const dt_limit_reason MAX_DT_NO_REASON_INFO = 1;

//Top level limits
const dt_limit_reason MAX_DT_CONFIGURED = 101;
const dt_limit_reason MAX_DT_COSMOLOGY_OUTPUT_REACHED = 102;
const dt_limit_reason MAX_DT_DATA_DUMP_REACHED = 103;
const dt_limit_reason MAX_DT_INITIAL_DT_LIMITED = 104;
const dt_limit_reason MAX_DT_STOP_TIME_REACHED = 105;
const dt_limit_reason MAX_DT_TIME_ACTION_REACHED = 106;

//Level limits
const dt_limit_reason MAX_DT_CONFIGURED_REFINED = 202;
const dt_limit_reason MAX_DT_UPPER_LEVEL_SYNC = 203;
const dt_limit_reason MAX_DT_UPPER_LEVEL_PRESYNC = 204;

//Grid limits
const dt_limit_reason MAX_DT_ACCELERATION_LIMITED = 301;
const dt_limit_reason MAX_DT_BURNING_ENERGY_LIMITED = 302;
const dt_limit_reason MAX_DT_CONDUCTION_LIMITED = 303;
const dt_limit_reason MAX_DT_COOLING_LIMITED = 304;
const dt_limit_reason MAX_DT_CR_LIMITED = 305;
const dt_limit_reason MAX_DT_EXPANSION_LIMITED = 306;
const dt_limit_reason MAX_DT_FREE_FALL_LIMITED = 307;
const dt_limit_reason MAX_DT_GASDRAG_LIMITED = 308;
const dt_limit_reason MAX_DT_HYDRO_LIMITED = 309;
const dt_limit_reason MAX_DT_HYDRO_RK_LIMITED = 310;
const dt_limit_reason MAX_DT_INTERNAL_ENERGY_GROWTH_LIMITED = 311;
const dt_limit_reason MAX_DT_MHD_LI_LIMITED = 312;
const dt_limit_reason MAX_DT_MHD_RK_LIMITED = 313;
const dt_limit_reason MAX_DT_RAD_PRESSURE_LIMITED = 314;
const dt_limit_reason MAX_DT_RAD_TRANSFER_LIMITED = 315;
const dt_limit_reason MAX_DT_SAFETY_VELOCITY_LIMITED = 316;
const dt_limit_reason MAX_DT_TOTAL_ENERGY_GROWTH_LIMITED = 317;
const dt_limit_reason MAX_DT_VISCOUS_LIMITED = 318;

//grid limits from particles
const dt_limit_reason MAX_DT_PARTICLES_LIMITED = 401;

struct DtLimitInfo
{
	float dt; //
	dt_limit_reason reason; //
	size_t dtIndex; // cell index in a grid context or another index
	int level; // hierarchy level
	int ProcessorNumber; // mpi process rank
	int GridID; // hierarchy grid id
	grid* Grid; //
	size_t ijk[MAX_DIMENSION]; //
	FLOAT xyz[MAX_DIMENSION]; //
	float uvw[MAX_DIMENSION]; //
	float rho, internalE;
};

#define MAX_DT_NEW_LOW_LIMIT tiny_number
#define MAX_DT_KEEP_SO_FAR ((dtNew <= MAX_DT_NEW_LOW_LIMIT) || (0 < *dtSoFar && *dtSoFar <= dtNew))

inline void LimitDt(float* const dtSoFar, const float dtNew)
{
	if(MAX_DT_KEEP_SO_FAR)
		return;

	*dtSoFar = dtNew;
}

inline void LimitDt(float* const dtSoFar, const float dtNew, size_t* dtIndexSoFar, const size_t dtIndexNew)
{
	if(MAX_DT_KEEP_SO_FAR)
		return;

	*dtSoFar = dtNew;
	*dtIndexSoFar = dtIndexNew;
}

void LimitDt(float* const dtSoFar, const float dtNew, DtLimitInfo* const infoSoFar, const DtLimitInfo* const infoNew);
void LimitDt(float* const dtSoFar, const float dtNew, DtLimitInfo* const infoSoFar, const dt_limit_reason reasonNew,
	const size_t dtIndexNew);

const char* const getReasonText(const dt_limit_reason reason);

int snlprintTimeStepLimitInfo(char* const s, const size_t size, size_t* const length, const DtLimitInfo* const info);

#endif /* LIMIT_TIME_STEP_H */
