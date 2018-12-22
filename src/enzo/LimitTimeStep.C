#include "myenzoutils.h"
#include "LimitTimeStep.h"
#include "global_data.h"
#include "Grid.h"

void LimitDt(float* const dtSoFar, const float dtNew, DtLimitInfo* const infoSoFar, const DtLimitInfo* const infoNew)
{
	if(MAX_DT_KEEP_SO_FAR)
		return;

	*dtSoFar = dtNew;
	if(infoSoFar && infoNew)
		*infoSoFar = *infoNew;
}

void LimitDt(float* const dtSoFar, const float dtNew, DtLimitInfo* const infoSoFar, const dt_limit_reason reasonNew,
	const size_t dtIndexNew)
{
	if(MAX_DT_KEEP_SO_FAR)
		return;

	*dtSoFar = dtNew;
	if(infoSoFar)
	{
		infoSoFar->dt = dtNew;
		infoSoFar->dtIndex = dtIndexNew;
		infoSoFar->reason = reasonNew;
	}
}

const char* const getReasonText(const dt_limit_reason reason)
{
	switch(reason)
	{
	case MAX_DT_ACCELERATION_LIMITED:
		return "MAX_DT_ACCELERATION_LIMITED";
	case MAX_DT_HYDRO_LIMITED:
		return "MAX_DT_HYDRO_LIMITED";
	case MAX_DT_HYDRO_RK_LIMITED:
		return "MAX_DT_HYDRO_RK_LIMITED";
	case MAX_DT_BURNING_ENERGY_LIMITED:
		return "MAX_DT_BURNING_ENERGY_LIMITED";
	case MAX_DT_CONDUCTION_LIMITED:
		return "MAX_DT_CONDUCTION_LIMITED";
	case MAX_DT_CONFIGURED_REFINED:
		return "CONFIGURED_LIMIT_REFINED";
	case MAX_DT_CONFIGURED:
		return "MAX_DT_CONFIGURED";
	case MAX_DT_COOLING_LIMITED:
		return "MAX_DT_COOLING_LIMITED";
	case MAX_DT_COSMOLOGY_OUTPUT_REACHED:
		return "MAX_DT_COSMOLOGY_OUTPUT_REACHED";
	case MAX_DT_CR_LIMITED:
		return "MAX_DT_CR_LIMITED";
	case MAX_DT_DATA_DUMP_REACHED:
		return "DATA_DUMP_REACHED";
	case MAX_DT_EXPANSION_LIMITED:
		return "MAX_DT_EXPANSION_LIMITED";
	case MAX_DT_FREE_FALL_LIMITED:
		return "MAX_DT_FREE_FALL_LIMITED";
	case MAX_DT_GASDRAG_LIMITED:
		return "MAX_DT_GASDRAG_LIMITED";
	case MAX_DT_INITIAL_DT_LIMITED:
		return "MAX_DT_INITIAL_DT_LIMITED";
	case MAX_DT_INTERNAL_ENERGY_GROWTH_LIMITED:
		return "MAX_DT_INTERNAL_ENERGY_GROWTH_LIMITED";
	case MAX_DT_NO_REASON_INFO:
		return "LIMITED_BY_ANOTHER_REASON";
	case MAX_DT_MHD_LI_LIMITED:
		return "MHD_LI_LIMITED";
	case MAX_DT_MHD_RK_LIMITED:
		return "MHD_RK_LIMITED";
	case MAX_DT_NULL_REASON:
		return "MAX_DT_NULL_REASON";
	case MAX_DT_PARTICLES_LIMITED:
		return "MAX_DT_PARTICLES_LIMITED";
	case MAX_DT_RAD_PRESSURE_LIMITED:
		return "MAX_DT_RAD_PRESSURE_LIMITED";
	case MAX_DT_RAD_TRANSFER_LIMITED:
		return "RAD_TRANSFER_LIMITED";
	case MAX_DT_STOP_TIME_REACHED:
		return "MAX_DT_STOP_TIME_REACHED";
	case MAX_DT_TIME_ACTION_REACHED:
		return "MAX_DT_TIME_ACTION_REACHED";
	case MAX_DT_TOTAL_ENERGY_GROWTH_LIMITED:
		return "MAX_DT_TOTAL_ENERGY_GROWTH_LIMITED";
	case MAX_DT_UPPER_LEVEL_PRESYNC:
		return "MAX_DT_UPPER_LEVEL_PRESYNC";
	case MAX_DT_UPPER_LEVEL_SYNC:
		return "MAX_DT_UPPER_LEVEL_SYNC";
	case MAX_DT_VISCOUS_LIMITED:
		return "MAX_DT_VISCOUS_LIMITED";
	default:
		return "MAX_DT_UNKNOWN_REASON";
	}
}

size_t sprintInfo(char* const s, const DtLimitInfo* const info)
{
	char* s2 = s;
	s2 += sprintf(s2, "%s", getReasonText(info->reason));
	if(info->reason == MAX_DT_NO_REASON_INFO)
	{
		//no details
	}
	else if(info->reason >= 200)
	{
		s2 += sprintf(s2, " level=%" ISYM, info->level);
		grid* g = info->Grid;
		if(g)
		{
			int n = g->GetGridRank();
			s2 += sprintfvec(s2, " cell=(", "%" ESYM, ",", ")", info->xyz, n, false, false);
			s2 += sprintfvec(s2, "[", "%" ISYM, ",", "]", info->ijk, n, false, false);
			s2 += sprintf(s2, "[%" ISYM "] grid[%" ISYM "]:", info->dtIndex, g->GetGridID());
			s2 += sprintfvec(s2, "(", "%" ESYM, ",", ")..", g->GetGridLeftEdge(), n, false, false);
			s2 += sprintfvec(s2, "(", "%" ESYM, ",", ")", g->GetGridRightEdge(), n, false, false);
			s2 += sprintfvec(s2, "/[", "%" ISYM, "x", "", g->GetGridDimension(), n, false, false);
			s2 += sprintf(s2, " = %lld]", g->GetGridSize());
		}
	}
	else if(info->reason >= 100)
	{
		s2 += sprintf(s2, " level=%" ISYM, info->level);
	}
	else
	{
		s2 += sprintf(s2, " index=%" ISYM, info->dtIndex);
	}

	return s2 - s;
}
