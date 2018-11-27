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

//in MHD_ObliqueRoutines.C
void RotateVector(float * Vector, float * Normal);
int SetupNormal(float Normal[], float MHDBlastCenter[3], TopGridData & MetaData);

int MHDProfileInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData,
		ExternalBoundary &Exterior)
{

	//
	//
	// Labels and Units.  (For IO.)
	//

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
	char *Density_56NiName = "Density_56Ni";   //[BH]
	char *Density_56NiUnits = NULL; //"g/cm**3";        //[BH]

	int i = 0, j = 0;

	DataLabel[i++] = DensName;
	DataUnits[j++] = NULL; //"g/cm**3";
	if (EquationOfState == 0)
	{
		DataLabel[i++] = TEName;
		DataUnits[j++] = NULL; //"erg/g";
	}
	DataLabel[i++] = Vel1Name;
	DataUnits[j++] = NULL; //"cm/s";
	DataLabel[i++] = Vel2Name;
	DataUnits[j++] = NULL; //"cm/s";
	DataLabel[i++] = Vel3Name;
	DataUnits[j++] = NULL; //"cm/s";
	if (UseMHD)
	{
		DataUnits[i] = NULL;
		DataLabel[i++] = BxName;
		DataUnits[i] = NULL;
		DataLabel[i++] = ByName;
		DataUnits[i] = NULL;
		DataLabel[i++] = BzName;
		DataUnits[i] = NULL;
		DataLabel[i++] = PhiName;
	}

	if (DualEnergyFormalism)
	{
		char *GEName = "GasEnergy";
		DataLabel[i++] = GEName;
		DataUnits[j++] = NULL; //"erg/g";
	}

	if (WritePotential)
	{
		DataLabel[i++] = GPotName;
		DataUnits[j++] = NULL;
	}

	if (UseBurning)
	{                         //[BH]
		DataUnits[i] = Density_56NiUnits;   //[BH]
		DataLabel[i++] = Density_56NiName;    //[BH]
//    DataUnits[i  ] = QInstantaneousUnits; //[BH]
//    DataLabel[i++] = QInstantaneousName;  //[BH]
//    DataUnits[i  ] = QCumulativeUnits;    //[BH]
//    DataLabel[i++] = QCumulativeName;     //[BH]
	}                                       //[BH]

	MHDLabel[0] = "BxF";
	MHDLabel[1] = "ByF";
	MHDLabel[2] = "BzF";

	MHDeLabel[0] = "Ex";
	MHDeLabel[1] = "Ey";
	MHDeLabel[2] = "Ez";

	MHDUnits[0] = "None";
	MHDUnits[1] = "None";
	MHDUnits[2] = "None";

	MHDeUnits[0] = "None";
	MHDeUnits[1] = "None";
	MHDeUnits[2] = "None";

	// General control variable
	int dim;

	// Parameters and their defaults.
	char line[MAX_LINE_LENGTH];
	int ret = 0, GasFlag = 0, Pflag = 0, TotalFlag = 0;
	int ObsFlag = 0;
	int RefineOnStartup = FALSE;

//  float DensityA = 1.0666;
//  float DensityB = 1.0;;
//  float GasEnergyA = 3.666;
//  float GasEnergyB = 1000.666;
//  float TotalEnergyA = 1.0;
//  float TotalEnergyB = 1.0;

//  float Rho_56Ni_InitialA = 0; //[BH]
//  float Rho_56Ni_InitialB = 1; //[BH]
	float InternalEnergy_InitialA = 0; //[BH] TODO: Perhaps come up with some initial values
	float InternalEnergy_InitialB = 0; //[BH]       for the intenal energy which make some sense.

	float MHDBlastSubgridLeft[3] =
	{ DomainLeftEdge[0], DomainLeftEdge[1], DomainLeftEdge[2] };
	float MHDBlastSubgridRight[3] =
	{ DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2] };

	//Obsolete variable names.
//  float Pressure0, Pressure1;
//  float B0[3],B1[3],Energy0, Energy1;
//  float Density0,Density1, GasEnergy0, GasEnergy1, TotalEnergy0,TotalEnergy1;

	char profileFileName[MAX_LINE_LENGTH];
	char profileFormat[16]; //"PLAIN" or "PAH1" or "PAH2"
	char profileType[16]; //"RADIAL"
	char radiusColumnName[32]; //"radius"
	char densityColumnName[32]; //"density"
	char temperatureColumnName[32]; //"temperature"
	float burningTemperature = 0;
	float burningRadius = 0;
	float profileAtTime = -1;

	*profileFileName = *profileFormat = *profileType = '\0';
	*radiusColumnName = *densityColumnName = *temperatureColumnName = '\0';

	//
	// Read Parameter File.
	//

	while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
	{

		ret = 0;
//		Use PSYM or FSYM for floats, ISYM for ints.
		ret += sscanf(line, "ProfileFileName = %s", &profileFileName);
		ret += sscanf(line, "ProfileType = %s", &profileType);
		ret += sscanf(line, "ProfileFormat = %s", &profileFormat);
		ret += sscanf(line, "RadiusColumnName = %s", &radiusColumnName);
		ret += sscanf(line, "DensityColumnName = %s", &densityColumnName);
		ret += sscanf(line, "TemperatureColumnName = %s", &temperatureColumnName);
		ret += sscanf(line, "BurningTemperature = %lf", &burningTemperature);
		ret += sscanf(line, "BurningRadius = %lf", &burningRadius);
		ret += sscanf(line, "ProfileAtTime = %lf", &profileAtTime);

//    ret += sscanf(line, "MHDBlastInitStyle = %"ISYM"", &InitStyle);
//    ret += sscanf(line, "Rho_56Ni_InitialA = %"FSYM"", &Rho_56Ni_InitialA);         //[BH]
//    ret += sscanf(line, "Rho_56Ni_InitialB = %"FSYM"", &Rho_56Ni_InitialB);         //[BH]
	ret += sscanf(line, "InternalEnergy_InitialA = %lf", &InternalEnergy_InitialA);//[BH]
	ret += sscanf(line, "InternalEnergy_InitialB = %lf", &InternalEnergy_InitialB);//[BH]
}
//line loop

//Re scale the subgrid edges to line up with the parent grid.
// nCellsL and nCellsR are the number of cells from the domain left edge.

int nCellsL[3], nCellsR[3];
int nCells[3] =
{ 0, 0, 0 };
for (dim = 0; dim < 3; dim++)
{
	nCellsL[dim] = nint(
			(MHDBlastSubgridLeft[dim] - DomainLeftEdge[dim]) / (DomainRightEdge[dim] - DomainLeftEdge[dim])
					* MetaData.TopGridDims[dim]);

	MHDBlastSubgridLeft[dim] = max(
			nCellsL[dim] * (DomainRightEdge[dim] - DomainLeftEdge[dim]) / MetaData.TopGridDims[dim],
			DomainLeftEdge[dim]);

	nCellsR[dim] = nint(
			(MHDBlastSubgridRight[dim] - DomainLeftEdge[dim]) / (DomainRightEdge[dim] - DomainLeftEdge[dim])
					* MetaData.TopGridDims[dim]);

	MHDBlastSubgridRight[dim] = min(
			nCellsR[dim] * (DomainRightEdge[dim] - DomainLeftEdge[dim]) / MetaData.TopGridDims[dim],
			DomainRightEdge[dim]);
	nCells[dim] = nint(
			(MHDBlastSubgridRight[dim] - MHDBlastSubgridLeft[dim]) / (DomainRightEdge[dim] - DomainLeftEdge[dim])
					* MetaData.TopGridDims[dim]);

}

if (RefineOnStartup == 1)
{
fprintf(stderr,"Subgrid Left %"GSYM" %"GSYM" %"GSYM"\n", MHDBlastSubgridLeft[0], MHDBlastSubgridLeft[1], MHDBlastSubgridLeft[2]);
fprintf(stderr,"Subgrid Right %"GSYM" %"GSYM" %"GSYM"\n",MHDBlastSubgridRight[0],MHDBlastSubgridRight[1],MHDBlastSubgridRight[2]);
fprintf(stderr,"nCells %"ISYM" %"ISYM" %"ISYM"\n", nCells[0], nCells[1], nCells[2]);
}

		// Long Dimension is used to convert the radius from Physical units to Grid Units;
		// We want the axis, though, so figure out which is the longest edge (in Grid Units)
		// then figure out which one it is.  A more elegant solution would be welcomed.

//  int LongDimension = 0;
//  LongDimension = (nCells[0] > nCells[1] ) ? nCells[0] : nCells[1];
//  LongDimension = (LongDimension > nCells[2] ) ? LongDimension : nCells[2];
//  for( dim=0; dim<3; dim++)
//    if( LongDimension == nCells[dim] ){
//      LongDimension = dim;
//      break;
//    }

//
// Calculate Total Energy.
//

if (Pflag > 0)
{
//      GasEnergyA = Pressure0/((Gamma-1)*DensityA);
//      GasEnergyB = Pressure1/((Gamma-1)*DensityB);
}

  //The variable stored is Gas+Kinetic+Magnetic Energy.
if (GasFlag > 0 || Pflag > 0)
{
//      Energy0 = GasEnergyA +
//	0.5*(VelocityA[0]*VelocityA[0] + VelocityA[1]*VelocityA[1] + VelocityA[2]*VelocityA[2])
//	+0.5*(BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2])/DensityA;
//      Energy1 = GasEnergyB +
//	0.5*(VelocityB[0]*VelocityB[0] + VelocityB[1]*VelocityB[1] + VelocityB[2]*VelocityB[2])
//	+0.5*(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2])/DensityB;
}

if (TotalFlag > 0)
{
//    Energy0=TotalEnergyA;
//    Energy1=TotalEnergyB;
}

  //
  // Initialize the top grid.  Cant' decide if I want a uniform grid here or MHDBlastInitialize.
  //
printf("MHDProfileInitialize ----------------------------------- 1\n");

if (TopGrid.GridData->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType, radiusColumnName,
	densityColumnName, temperatureColumnName, burningTemperature, burningRadius, profileAtTime) == FAIL)
ENZO_FAIL("MHDProfileInitialize:  Error in MHDProfileInitializeGrid.");

printf("MHDProfileInitialize ----------------------------------- 2\n");
  //	//
  // MagneticField Initialize.
  //
  // Generate Hierarchy.
  //
if (RefineOnStartup == 1)
{
//Create as many subgrids as there are refinement levels
//needed to resolve the initial explosion region upon the start-up.

HierarchyEntry ** Subgrid;
if (MaximumRefinementLevel > 0)
	Subgrid = new HierarchyEntry*[MaximumRefinementLevel];

//
//Create new HierarchyEntries.  Note that 'lev' loops are only for the SUBGRIDS.
//

int lev;
int NumberOfSubgridZones[3], SubgridDims[3];

for (lev = 0; lev < MaximumRefinementLevel; lev++)
	Subgrid[lev] = new HierarchyEntry;

for (lev = 0; lev < MaximumRefinementLevel; lev++)
{
	//Calculate number of cells on this level.

	for (dim = 0; dim < MetaData.TopGridRank; dim++)
		NumberOfSubgridZones[dim] = nCells[dim] * POW(RefineBy, lev + 1);

	fprintf(stderr,"uncle MHDBlast:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
			NumberOfSubgridZones[0]);

	if (NumberOfSubgridZones[0] > 0)
	{
		// fill them out

		if (lev == 0)
			TopGrid.NextGridNextLevel = Subgrid[0];
		Subgrid[lev]->NextGridThisLevel = NULL;
		if (lev == MaximumRefinementLevel - 1)
			Subgrid[lev]->NextGridNextLevel = NULL;
		else
			Subgrid[lev]->NextGridNextLevel = Subgrid[lev + 1];
		if (lev == 0)
			Subgrid[lev]->ParentGrid = &TopGrid;
		else
			Subgrid[lev]->ParentGrid = Subgrid[lev - 1];

		//  compute the dimensions and left/right edges for the subgrid

		for (dim = 0; dim < MetaData.TopGridRank; dim++)
		{
			SubgridDims[dim] = NumberOfSubgridZones[dim] + 2 * NumberOfGhostZones;
		}

		// create a new subgrid and initialize it

		Subgrid[lev]->GridData = new grid;
		Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
		Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims, MHDBlastSubgridLeft,
				MHDBlastSubgridRight, 0);

		if (Subgrid[lev]->GridData->MHDProfileInitializeGrid(profileFileName, profileFormat, profileType,
				radiusColumnName, densityColumnName, temperatureColumnName, burningTemperature, burningRadius,
				profileAtTime) == FAIL)
			ENZO_FAIL("MHDProfileInitialize: Error in MHDProfileInitializeGrid.");

	}    //NumberOfSubgridZones > 0
	else
	{
		printf("single grid start-up.\n");
	}

}    //level

// Make sure each grid has the best data with respect to the finer grids.
// This projection juggle is to ensure that, regardless of how the hierarchy is evolved, the field gets projected
// properly here.

int MHD_ProjectEtmp = MHD_ProjectE;
int MHD_ProjectBtmp = MHD_ProjectB;
MHD_ProjectE = FALSE;
MHD_ProjectB = TRUE;

for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
	printf("projecting ----------> %d\n", lev);;
	if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(*(Subgrid[lev - 1]->GridData)) == FAIL)
		ENZO_FAIL("Error in ProjectSolutionToParentGrid.");

// set up the root grid

if (MaximumRefinementLevel > 0)
{
	if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData)) == FAIL)
		ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
}

// Put the projection options back to the inital.
MHD_ProjectE = MHD_ProjectEtmp;
MHD_ProjectB = MHD_ProjectBtmp;

}    //RefineOnStartup

return SUCCESS;
}

