/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include "ParameterControl/ParameterControl.h"
extern Configuration Param;

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

/* Set default parameter values. */

const char config_collapse_test_defaults[] = 
"### COLLAPSE TEST INITIALIZATION DEFAULTS ###\n"
"\n"
"Problem: {\n"
"    CollapseTest: {\n"
"        RefineAtStart = True;\n"
"        UseParticles  = False;\n"
"        ParticleMeanDensity = 0.0;\n"
"        UseColour = False;\n"
"        UseMetals = False;\n"
"        InitialTemperature = 1000;\n"
"        InitialDensity     = 1.0;\n"
"\n"
"        Spheres = [\"Sphere1\"];\n"
"\n"
"        Sphere1: {\n"
"            Density = 1.0;\n"
"            Temperature = 1.0;\n"
"            Velocity = [0.0, 0.0, 0.0];\n"
"            UniformVelocity = 0.0;\n"
"            FracKeplerianRot = 0.0;\n"
"            Turbulence = 0.0;\n"
"            Dispersion = 0.0;\n"
"            CutOff = 6.5;\n"
"            Ang1 = 0.0;\n"
"            Ang2 = 0.0;\n"
"            Metallicity = 0.0;\n"
"            NumShells = 1;\n"
"            InitialLevel = 0;\n"
"            Type = 0;\n"
"            Radius = 1.0;\n"
"            CoreRadius = 0.1;\n"
"            Position = [-99999.0, -99999.0, -99999.0];\n"
"        };\n"
"    };\n"
"};\n";
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int CollapseTestInitialize(FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "colour";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *MetalName = "Metal_Density";

  /* declarations */

  int   dim, level, sphere, i;

  /* set default parameters */

  int CollapseTestNumberOfSpheres = 1;
  int CollapseTestRefineAtStart   = TRUE;
  int CollapseTestUseParticles    = FALSE;
  float CollapseTestParticleMeanDensity = FLOAT_UNDEFINED;
  int CollapseTestUseColour       = FALSE;
  int CollapseTestUseMetals       = FALSE;
  float CollapseTestInitialTemperature = 1000;
  float CollapseTestInitialDensity     = 1.0;
  float CollapseTestSphereDensity[MAX_SPHERES],
    CollapseTestSphereTemperature[MAX_SPHERES],
    CollapseTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
    CollapseTestUniformVelocity[MAX_DIMENSION],
    CollapseTestFracKeplerianRot[MAX_SPHERES],
    CollapseTestSphereTurbulence[MAX_SPHERES],
    CollapseTestSphereDispersion[MAX_SPHERES],
    CollapseTestSphereCutOff[MAX_SPHERES],
    CollapseTestSphereAng1[MAX_SPHERES],
    CollapseTestSphereAng2[MAX_SPHERES],
    CollapseTestSphereMetallicity[MAX_SPHERES];
  int CollapseTestSphereNumShells[MAX_SPHERES],
    CollapseTestSphereInitialLevel[MAX_SPHERES],
    CollapseTestSphereType[MAX_SPHERES];
  FLOAT CollapseTestSphereRadius[MAX_SPHERES],
    CollapseTestSphereCoreRadius[MAX_SPHERES],
    CollapseTestSpherePosition[MAX_SPHERES][MAX_DIMENSION];


  // This is how it should look eventually.
  //Param.UpdateDefaults(config_collapse_test_defaults);

  // We need to update the Sphere1.Position parameter in the defaults
  // settings, since it depends on runtime parameters (DomainLeftEdge,
  // DomainRightEdge). Super kludgy...

  FLOAT DefaultSphere1Position[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
      DefaultSphere1Position[dim] = 0.5*(DomainLeftEdge[dim] +
					 DomainRightEdge[dim]);
  //  Param.UpdateDefaults("Problem.CollapseTest.Sphere1.Position",DefaultSphere1Position);


  /* read input from file */

  Param.GetScalar(CollapseTestRefineAtStart, "Problem.CollapseTest.RefineAtStart");
  Param.GetScalar(CollapseTestUseParticles, "Problem.CollapseTest.UseParticles");
  Param.GetScalar(CollapseTestParticleMeanDensity, "Problem.CollapseTest.ParticleMeanDensity");
  Param.GetScalar(CollapseTestUseColour, "Problem.CollapseTest.UseColour");
  Param.GetScalar(CollapseTestUseMetals, "Problem.CollapseTest.UseMetals");
  Param.GetScalar(CollapseTestInitialTemperature, "Problem.CollapseTest.InitialTemperature");
  Param.GetScalar(CollapseTestInitialDensity, "Problem.CollapseTest.InitialDensity");
  Param.GetArray(CollapseTestUniformVelocity, "Problem.CollapseTest.UniformVelocity");
  
  int NumberOfSpheres = Param.Size("Problem.CollapseTest.Spheres");
  if (NumberOfSpheres > MAX_SPHERES-1) {
    ENZO_VFAIL("You've exceeded the maximum number of CollapseTest spheres (%d)!\n",MAX_SPHERES)
      }
  CollapseTestNumberOfSpheres = NumberOfSpheres;
  
  char SphereNames[MAX_LINE_LENGTH][MAX_SPHERES];
  Param.GetArray(SphereNames,"Problem.CollapseTest.Spheres");
  
  for (i = 0; i < NumberOfSpheres; i++) {
    Param.GetScalar(CollapseTestSphereType[i], "Problem.CollapseTest.%s.Type",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereRadius[i], "Problem.CollapseTest.%s.Radius",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereCoreRadius[i], "Problem.CollapseTest.%s.CoreRadius",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereDensity[i], "Problem.CollapseTest.%s.Density",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereTemperature[i], "Problem.CollapseTest.%s.Temperature",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereMetallicity[i], "Problem.CollapseTest.%s.Metallicity",SphereNames[i]);
    Param.GetArray(CollapseTestSpherePosition[i], "Problem.CollapseTest.%s.Position",SphereNames[i]);
    Param.GetArray(CollapseTestSphereVelocity[i], "Problem.CollapseTest.%s.Velocity",SphereNames[i]);
    Param.GetScalar(CollapseTestFracKeplerianRot[i], "Problem.CollapseTest.%s.FracKeplerianRot",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereTurbulence[i], "Problem.CollapseTest.%s.Turbulence",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereDispersion[i], "Problem.CollapseTest.%s.Dispersion",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereCutOff[i], "Problem.CollapseTest.%s.CutOff",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereAng1[i], "Problem.CollapseTest.%s.Ang1",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereAng2[i], "Problem.CollapseTest.%s.Ang2",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereNumShells[i], "Problem.CollapseTest.%s.NumShells",SphereNames[i]);
    Param.GetScalar(CollapseTestSphereInitialLevel[i], "Problem.CollapseTest.%s.InitialLevel",SphereNames[i]);
  } 

  
  /* set up grid */

  if (TopGrid.GridData->CollapseTestInitializeGrid(
	     CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
	     CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
	     CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
	     CollapseTestSpherePosition, CollapseTestSphereVelocity,
             CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
	     CollapseTestSphereDispersion,
             CollapseTestSphereCutOff, CollapseTestSphereAng1,
             CollapseTestSphereAng2, CollapseTestSphereNumShells,
	     CollapseTestSphereType, CollapseTestUseParticles,
	     CollapseTestParticleMeanDensity,
             CollapseTestUniformVelocity, CollapseTestUseColour,
	     CollapseTestUseMetals,
             CollapseTestInitialTemperature, CollapseTestInitialDensity,
	     0) == FAIL) {
    ENZO_FAIL("Error in CollapseTestInitializeGrid.");
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested and there are no manual settings of the refinement
     of spheres, refine the grid to the desired level. */

  int MaxInitialLevel = 0;
  for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++)
    MaxInitialLevel = max(MaxInitialLevel, CollapseTestSphereInitialLevel[sphere]);

  if (CollapseTestRefineAtStart) {

    /* If the user specified an initial refinement level for a sphere,
       then manually create the hierarchy first. */

    if (MaxInitialLevel > 0) {

      int lev, max_level;
      float dx;
      HierarchyEntry **Subgrid;
      int NumberOfSubgridDims[MAX_DIMENSION];
      FLOAT ThisLeftEdge[MAX_DIMENSION], ThisRightEdge[MAX_DIMENSION];

      for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
	
	max_level = CollapseTestSphereInitialLevel[sphere];
	if (max_level > 0) {

	  Subgrid = new HierarchyEntry*[max_level];
	  for (lev = 0; lev < max_level; lev++)
	    Subgrid[lev] = new HierarchyEntry;

	  for (lev = 0; lev < max_level; lev++) {
	    
	    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	      dx = 1.0 / float(MetaData.TopGridDims[dim]) / POW(RefineBy, lev);
	      ThisLeftEdge[dim] = CollapseTestSpherePosition[sphere][dim] -
		0.5 * CollapseTestSphereRadius[sphere] - 2*dx;  // plus some buffer
	      ThisLeftEdge[dim] = nint(ThisLeftEdge[dim] / dx) * dx;
	      ThisRightEdge[dim] = CollapseTestSpherePosition[sphere][dim] +
		0.5 * CollapseTestSphereRadius[sphere] + 2*dx;
	      ThisRightEdge[dim] = nint(ThisRightEdge[dim] / dx) * dx;
	      NumberOfSubgridDims[dim] = 
		nint((ThisRightEdge[dim] - ThisLeftEdge[dim]) / 
		     (DomainRightEdge[dim] - DomainLeftEdge[dim]) / dx);		
	    } // ENDFOR dims

	    if (debug)
	      printf("CollapseTest:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n",
		     lev+1, NumberOfSubgridDims[0]);
	    
	    if (NumberOfSubgridDims[0] > 0) {

	      // Insert into AMR hierarchy
	      if (lev == 0) {
		Subgrid[lev]->NextGridThisLevel = TopGrid.NextGridNextLevel;
		TopGrid.NextGridNextLevel = Subgrid[lev];
		Subgrid[lev]->ParentGrid = &TopGrid;
	      } else {
		Subgrid[lev]->NextGridThisLevel = NULL;
		Subgrid[lev]->ParentGrid = Subgrid[lev-1];
	      }
	      if (lev == max_level-1)
		Subgrid[lev]->NextGridNextLevel = NULL;
	      else
		Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];

	      // Create grid
	      for (dim = 0; dim < MetaData.TopGridRank; dim++)
		NumberOfSubgridDims[dim] += 2*DEFAULT_GHOST_ZONES;
	      Subgrid[lev]->GridData = new grid;
	      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
	      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, 
						  NumberOfSubgridDims,
						  ThisLeftEdge,
						  ThisRightEdge, 0);


	      if (Subgrid[lev]->GridData->CollapseTestInitializeGrid(
	          CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
		  CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
		  CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
		  CollapseTestSpherePosition, CollapseTestSphereVelocity,
		  CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
		  CollapseTestSphereDispersion,
		  CollapseTestSphereCutOff, CollapseTestSphereAng1,
		  CollapseTestSphereAng2, CollapseTestSphereNumShells,
		  CollapseTestSphereType, CollapseTestUseParticles,
		  CollapseTestParticleMeanDensity,
		  CollapseTestUniformVelocity, CollapseTestUseColour,
		  CollapseTestUseMetals,
		  CollapseTestInitialTemperature, CollapseTestInitialDensity,
		  lev-1) == FAIL) {
		ENZO_FAIL("Error in CollapseTestInitializeGrid.");
	      }
	      
	    } // ENDIF zones exist
	  } // ENDFOR levels
	} // ENDIF max_level > 0
      } // ENDFOR spheres
    } // ENDIF MaxInitialLevel > 0

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    if (MaxInitialLevel == 0) {
      for (level = 0; level < MaximumRefinementLevel; level++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	  ENZO_FAIL("Error in RebuildHierarchy.");
	}
	if (LevelArray[level+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level+1];
	while (Temp != NULL) {
	  if (Temp->GridData->CollapseTestInitializeGrid(
		 CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
		 CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
		 CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
		 CollapseTestSpherePosition, CollapseTestSphereVelocity,
		 CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
		 CollapseTestSphereDispersion,
		 CollapseTestSphereCutOff, CollapseTestSphereAng1,
		 CollapseTestSphereAng2, CollapseTestSphereNumShells,
		 CollapseTestSphereType, CollapseTestUseParticles,
		 CollapseTestParticleMeanDensity,
		 CollapseTestUniformVelocity, CollapseTestUseColour,
		 CollapseTestUseMetals,
		 CollapseTestInitialTemperature, CollapseTestInitialDensity,
		 level+1) == FAIL) {
	    ENZO_FAIL("Error in CollapseTestInitializeGrid.");
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      } // end: loop over levels
    } // ENDELSE manually set refinement levels

      /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
			      *LevelArray[level-1]->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (CollapseTestRefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = (char*) GEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
  if (MultiSpecies) {
    DataLabel[count++] = (char*) ElectronName;
    DataLabel[count++] = (char*) HIName;
    DataLabel[count++] = (char*) HIIName;
    DataLabel[count++] = (char*) HeIName;
    DataLabel[count++] = (char*) HeIIName;
    DataLabel[count++] = (char*) HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = (char*) HMName;
      DataLabel[count++] = (char*) H2IName;
      DataLabel[count++] = (char*) H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = (char*) DIName;
      DataLabel[count++] = (char*) DIIName;
      DataLabel[count++] = (char*) HDIName;
    }
  }  // if Multispecies
  if (CollapseTestUseColour)
    DataLabel[count++] = (char*) ColourName;
  if (CollapseTestUseMetals)
    DataLabel[count++] = (char*) MetalName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CollapseTestNumberOfSpheres    = %"ISYM"\n",
	    CollapseTestNumberOfSpheres);
    fprintf(Outfptr, "CollapseTestRefineAtStart      = %"ISYM"\n",
	    CollapseTestRefineAtStart);
    fprintf(Outfptr, "CollapseTestUseParticles       = %"ISYM"\n",
	    CollapseTestUseParticles);
    fprintf(Outfptr, "CollapseTestUseColour          = %"ISYM"\n",
	    CollapseTestUseColour);
    fprintf(Outfptr, "CollapseTestUseMetals          = %"ISYM"\n",
	    CollapseTestUseMetals);
    fprintf(Outfptr, "CollapseTestInitialTemperature = %"FSYM"\n",
	    CollapseTestInitialTemperature);
    fprintf(Outfptr, "CollapseTestInitialDensity     = %"FSYM"\n",
	    CollapseTestInitialDensity);
    fprintf(Outfptr, "CollapseTestUniformVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
	    CollapseTestUniformVelocity[0], CollapseTestUniformVelocity[1],
	    CollapseTestUniformVelocity[2]);
    for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
      fprintf(Outfptr, "CollapseTestSphereType[%"ISYM"] = %"ISYM"\n", sphere,
	      CollapseTestSphereType[sphere]);
      fprintf(Outfptr, "CollapseTestSphereRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereCoreRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereCoreRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereDensity[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereDensity[sphere]);
      fprintf(Outfptr, "CollapseTestSphereTemperature[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereTemperature[sphere]);
      fprintf(Outfptr, "CollapseTestSphereMetallicity[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereMetallicity[sphere]);
      fprintf(Outfptr, "CollapseTestSpherePosition[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSpherePosition[sphere]);
      fprintf(Outfptr, "CollapseTestSphereVelocity[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSphereVelocity[sphere]);
      fprintf(Outfptr, "CollapseTestFracKeplerianRot[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestFracKeplerianRot[sphere]);
      fprintf(Outfptr, "CollapseTestSphereTurbulence[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereTurbulence[sphere]);
      fprintf(Outfptr, "CollapseTestSphereCutOff[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereCutOff[sphere]);
      fprintf(Outfptr, "CollapseTestSphereAng1[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereAng1[sphere]);
      fprintf(Outfptr, "CollapseTestSphereAng2[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereAng2[sphere]);
      fprintf(Outfptr, "CollapseTestSphereNumShells[%"ISYM"] = %"ISYM"\n", sphere,
              CollapseTestSphereNumShells[sphere]);
    }
  }

  return SUCCESS;

}
