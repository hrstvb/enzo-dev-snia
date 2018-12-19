/***********************************************************************
 /
 /  SET THE TIMESTEP FOR A LEVEL
 /
 /  written by: Greg Bryan
 /  date:       November, 1994
 /  modified1:  Matthew Turk, split off
 /  date:       June 2009
 /
 /  PURPOSE:
 /       Determine the timestep for this iteration of the loop.
 /
 ************************************************************************/

#include "performance.h"
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
#include "LevelHierarchy.h"
#include "LimitTimeStep.h"

float CommunicationMinValue(float Value);

int SetLevelTimeStep(HierarchyEntry *Grids[], int NumberOfGrids, int level,
float *dtThisLevelSoFar, float *dtThisLevel,
float dtLevelAbove)
{
	float dt = huge_number;
	float dtGrid = huge_number;
	float dtConduction = huge_number;
	float dtRefined = huge_number;
	float dtSyncUpper = huge_number;
	float dtPreSyncUpper = huge_number;
	float dtGridIndex = MAX_DT_NO_INDEX_INFO;
	float dtConductionIndex = MAX_DT_NO_INDEX_INFO;
	float dtGrids;
	int grid1;
	struct DtLimitInfo gridsInfo;
	struct DtLimitInfo gridInfoNew;
	struct DtLimitInfo levelInfo;

	LCAPERF_START("SetLevelTimeStep"); // SetTimeStep()

	if(level == 0)
	{
		/* For root level, use dtLevelAbove. */
		*dtThisLevelSoFar = dt = dtLevelAbove;
		// This info is printed in EvolveHierarchy.
//		if(MyProcessorNumber == ROOT_PROCESSOR)
//			printf("Level[%" ISYM "](#%lld): dt = %" GSYM "\n", level, MyProcessorNumber, dt);
	}
	else
	{
		int my_isotropic_conduction = IsotropicConduction;
		int my_anisotropic_conduction = AnisotropicConduction;
		int dynamic_hierarchy_rebuild = ConductionDynamicRebuildHierarchy && level >= ConductionDynamicRebuildMinLevel
				&& dtRebuildHierarchy[level] <= 0.0;

		if(dynamic_hierarchy_rebuild)
		{
			/* Calculate timestep without conduction and get conduction separately later. */
			IsotropicConduction = AnisotropicConduction = FALSE;
		}

		/* Compute the mininum timestep for all grids. */

		for(grid1 = 0; grid1 < NumberOfGrids; grid1++)
		{
			dtGrid = Grids[grid1]->GridData->ComputeTimeStep(&gridInfoNew);
			//*dtThisLevel = min(*dtThisLevel, dtGrid);
			LimitDt(&dt, dtGrid, &gridsInfo, &gridInfoNew);
		}
		dt = CommunicationMinValue(dt);

		/* Compute conduction timestep and use to set the number
		 of iterations without rebuiding the hierarchy. */

		if(dynamic_hierarchy_rebuild)
		{
			if(my_isotropic_conduction || my_anisotropic_conduction)
			{
				/* Return conduction parameters to original values. */
				IsotropicConduction = my_isotropic_conduction;
				AnisotropicConduction = my_anisotropic_conduction;

				float dt_cond_temp;
				for(grid1 = 0; grid1 < NumberOfGrids; grid1++)
				{
					if(Grids[grid1]->GridData->ComputeConductionTimeStep(dt_cond_temp) == FAIL)
						ENZO_FAIL("Error in ComputeConductionTimeStep.\n");
					dtConduction = min(dtConduction, dt_cond_temp);
				}
				dtConduction = CommunicationMinValue(dtConduction);
				dtConduction *= float(NumberOfGhostZones);  // for subcycling

				int my_cycle_skip = max(1, (int) (dt / dtConduction));
				dtRebuildHierarchy[level] = *dtThisLevel;
				if(debug)
					fprintf(stderr,
							"Conduction dt[%"ISYM"] = %"GSYM", will rebuild hierarchy in about %"ISYM" cycles.\n",
							level, dtConduction, my_cycle_skip);

				/* Set actual timestep correctly. */
				// *dtThisLevel = min(*dtThisLevel, dtConduction);
				LimitDt(&dt, dtConduction, &gridsInfo, MAX_DT_CONDUCTION_LIMITED, dtConductionIndex);
			}
			else
			{
				if(debug)
				{
					fprintf(stderr,
							"Isotrpic and anisotropic conductions are off. Dynamic hieararchy rebuid has no effect.\n");
				}
			}
		}

		dtGrids = dt; // Reduceall(min) result
		levelInfo = gridsInfo; // local values

#ifdef USE_DT_LIMIT
		//dtRefined = LevelZeroDeltaT/(4.0)/POW(RefineBy,level);

		//The following wouldn't result in the level0t for level = 0;
		dtRefined = /* MAX_TOP_LEVEL_DT??? */0.5 / (4.0) / POW(2.0, level);
		//Also this can be calculated from the upper level time step, i.g.
		//    dtRefined = dtLevelAboveLastTimeStep / RefineBy;
		//which is a more restrictive value since
		//    dtLebelAboe <= level0dt / RefineBy**level.

		//if ( dtActual < dtRefined )
		//  *dtThisLevel = dtRefined;
		LimitDt(&dt, dtRefined, &levelInfo, MAX_DT_CONFIGURED_MAX_DT_REFINED, MAX_DT_INDEX_NOT_APPLICABLE);
#endif

		/* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
		//if(*dtThisLevelSoFar + *dtThisLevel * 1.05 >= dtLevelAbove)
		//{
		//	*dtThisLevel = dtLevelAbove - *dtThisLevelSoFar;
		//	*dtThisLevelSoFar = dtLevelAbove;
		//}
		//else
		//	*dtThisLevelSoFar += *dtThisLevel;
		// The grids on this level must not time-step ahead of the
		// grids on the upper level, which may require to truncate
		// the calculated dt.  Even without truncation, we make
		// sure that the remaining time until syncing with the
		// upper level is not arbitrarily small, which will cause
		// the next time step to be truncated to an arbitrarily
		// small number.
		dtSyncUpper = dtLevelAbove - *dtThisLevelSoFar;
		if(dtSyncUpper <= 1.25 * dt)
		{
			if(dtSyncUpper > dt)
			{
				dtPreSyncUpper = 0.75 * dtSyncUpper;
				*dtThisLevelSoFar += dtPreSyncUpper;
				LimitDt(&dt, dtPreSyncUpper, &levelInfo, MAX_DT_UPPER_LEVEL_PRESYNC, MAX_DT_INDEX_NOT_APPLICABLE);
			}
			else
			{
				*dtThisLevelSoFar = dtLevelAbove;
				LimitDt(&dt, dtSyncUpper, &levelInfo, MAX_DT_UPPER_LEVEL_SYNC, MAX_DT_INDEX_NOT_APPLICABLE);
			}
		}
		else
		{
			*dtThisLevelSoFar += dt;
		}

		// Here dtGrids holds the minimum dt reduced from all mpi
		// ranks.  gridsInfo has the local (process's) values,
		// including the local dt.  If dt==dtGrids, then this
		// process has the grid from where the min dt originates.
		// In this case we print debug info from the current
		// process, instead the root process.  As a side effect
		// multiple printouts at the same level may appear if
		// another process happens to have a grid, which yielded
		// the same dt.  An alternative is to communicate the info
		// structure to the root process and print it from there.
		if(dtGrids == gridsInfo.dt)
		{
			// If dt had changed from dtGrids, it means that stricter
			// level limits had been applied.  In this case we print
			// level info, otherwise only the grid info.
			if(dtGrids != dt)
			{
				fprintf(stderr, "Level[%" ISYM "](#%lld): dt = %" GSYM "    dtTillSync = %" GSYM
				"    reason = %s    dtGrids = %" GSYM "    gridReason = %s\n",
						level, gridsInfo.ProcessorNumber, dt, dtSyncUpper, getReasonText(levelInfo.reason), dtGrids,
						getReasonText(gridsInfo.reason));
			}
			else
			{
				fprintf(stderr, "Level[%" ISYM "](#%lld): dt = %" GSYM "    dtTillSync = %" GSYM
				"    reason = %s\n",
						level, gridsInfo.ProcessorNumber, dt, dtSyncUpper, getReasonText(gridsInfo.reason));
			}
		}
	} // end if level==0

	//Set the output value.
	*dtThisLevel = dt;

	/* Set all grid's timestep to this minimum dt. */
	for(grid1 = 0; grid1 < NumberOfGrids; grid1++)
		Grids[grid1]->GridData->SetTimeStep(dt);

	if(debug)
		printf("Level[%"ISYM"]: dt = %" GSYM "  %" GSYM " (%" GSYM "/%" GSYM ")\n", level, *dtThisLevel, dtGrids,
				*dtThisLevelSoFar, dtLevelAbove);

	LCAPERF_STOP("SetLevelTimeStep"); // SetTimeStep()
	return SUCCESS;
}
