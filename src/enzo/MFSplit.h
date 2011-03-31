/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Frequency, Multi-species, Split Problem Class
/
/  written by: Daniel Reynolds
/  date:       June 2009
/  modified1:  
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit multi-frequency solve.
/
************************************************************************/

#ifdef TRANSFER
#ifndef MF_SPLIT_PROBLEM_DEFINED__
#define MF_SPLIT_PROBLEM_DEFINED__

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "InexactNewton.h"
#include "ImplicitProblemABC.h"
#include "FSProb.h"


class MFSplit : public virtual ImplicitProblemABC {

 private:
  
  // overall time spent in solver
  float RTtime;
  float timers[30];
  
  // Free-streaming problem solver, old/new time steps
  FSProb *FSSolve;
  float *Efree;

  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size (per frequency)
#ifdef USE_HYPRE
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
#endif

  // HYPRE Solver-specific data
  float  sol_tolerance;          // desired linear solver tolerance
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_rlxtype;            // relaxation type:
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  Eint32 SolvIndices[3][2];      // L/R edge indices of subdomain in global mesh
                                 // Note: these INCLUDE Dirichlet zones, even 
                                 //   though those are not included as active 
                                 //   data in the vectors or physics routines.
  int SolvOff[3];                // offset between HYPRE mesh and active mesh; 
                                 //   typically 0, but will be 1 for inclusion 
                                 //   of Dirichlet zones in HYPRE grid.

  // HYPRE interface temporary data
#ifdef USE_HYPRE
  HYPRE_StructMatrix P;          // holds radiation matrices
  HYPRE_StructVector rhsvec;     // holds radiation rhs vectors
  HYPRE_StructVector solvec;     // holds radiation solution vectors
#endif
  Eflt64 *matentries;            // holds matrix entries
  Eflt64 *rhsentries;            // holds rhs entries
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face (includes Dirichlet bdry zones)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dimension
  float EdgeVals[3][2];    // L/R edges of this proc's subdomain
  float *EBdryVals[3][2];  // boundary values for radiation BCs

  // time-stepping related data
  int initial_guess;   // parameter for setting the initial guess
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation time step size
  float mindt;         // minimum radiation time step size
  float dtfac[3];      // desired relative change in fields per step
  float dtnorm;        // norm choice for computing relative change:
                       //    0 -> max pointwise norm (default)
                       //   >0 -> rms p-norm over entire domain
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dtchem;        // chemistry + gas energy time step size
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  int LimType;         // flux limiter formulation:
                       //    0 -> standard Levermore-Pomraning limiter (LP, 1981)
                       //    1 -> rational approx. to LP limiter (LP, 1981)
                       //    2 -> Reynolds approx. to LP limiter
                       //    3 -> no limiter
                       //    4 -> ZEUS limiter (like 1, but no 'albedo')
  EnzoVector *sol;     // solution vector
  EnzoVector *U0;      // old time-level state
  EnzoVector *extsrc;  // temporary vector holding external forcing sources
  
  // problem defining data
  int iE1;             // EnzoVector index of first radiation frequency field
  int iE2;             // EnzoVector index of first radiation frequency field
  int iE3;             // EnzoVector index of first radiation frequency field
  int iec;             // EnzoVector index of first radiation frequency field
  int iHI;             // EnzoVector index of first radiation frequency field
  int iHeI;            // EnzoVector index of first radiation frequency field
  int iHeII;           // EnzoVector index of first radiation frequency field
  int Nchem;           // number of chemical species (non-negative integer)
  int Model;           // model choice, 0=>decoupled ODE test case
                       //               1=>case B HII recomb, no emissivity
                       //               4=>isothermal ionization, pt. src. rad.
  float NGammaDot;     // ionization strength (photons/sec)
  float EtaRadius;     // ionization source radius
  float EtaCenter[3];  // ionization source location

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  float rScale;        // FS radiation energy density scaling factor
  float eScale;        // specific energy correction scaling factor
  float nScale;        // species density scaling factor
  float fsUnits;       // free-streaming radiation unit conversion factor
  float rUnits;        // radiation energy density unit conversion factor
  float E1Units;       // nu0_HI monochromatic radiation scaling factor
  float E2Units;       // nu0_HeI monochromatic radiation scaling factor
  float E3Units;       // nu0_HeII monochromatic radiation scaling factor
  float eUnits;        // specific energy correction unit conversion factor
  float nUnits;        // species density unit conversion factor
  float rUnits0;       // radiation energy density unit conversion factor (old time)
  float E1Units0;      // nu0_HI monochromatic radiation scaling factor (old time)
  float E2Units0;      // nu0_HeI monochromatic radiation scaling factor (old time)
  float E3Units0;      // nu0_HeII monochromatic radiation scaling factor (old time)
  float nUnits0;       // species density unit conversion factor (old time)
  float DenUnits;      // density scaling factor
  float LenUnits;      // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor
  float DenUnits0;     // density scaling factor (old time)
  float LenUnits0;     // length scaling factor (old time)
  float chibar;        // int_nu0^infty chi(nu) d nu
  float chinuint;      // int_0^infty chi(nu) nu0/nu d nu

  // chemistry constants
  float HFrac;         // Fraction of matter composed of Hydrogen

  // storage for integrals over radiation spectrum
  float hnu0_HI;       // HI ionization threshold (eV)
  float hnu0_HeI;      // HeI ionization threshold (eV)
  float hnu0_HeII;     // HeII ionization threshold (eV)
  int ESpectrum;       // radiation spectrum type 
                       // (see MFSplit_ComputeRadiationIntegrals.src90)
  float *piHI;         // HI photo-ionization field
  float *piHeI;        // HeI photo-ionization field
  float *piHeII;       // HeII photo-ionization field
  float *GHI;          // HI photo-heating coefficient field
  float *GHeI;         // HeI photo-heating coefficient field
  float *GHeII;        // HeII photo-heating coefficient field

  // access to Enzo data
  float *vx;           // x0-directional velocity
  float *vy;           // x1-directional velocity
  float *vz;           // x2-directional velocity
  float *rho;          // density
  float *eh;           // fluid energy (total or internal)

  // stored arrays for increased efficiency
  float *eCorr;        // gas energy correction
  
  // private computation routines
  int EnforceBoundary(EnzoVector *vec);
  int SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, int freq, 
		  EnzoVector *u, EnzoVector *src, float dt);
  int ComputeRadiationIntegrals(EnzoVector *u);
  int EnforceRadiationBounds(float *Ef, float *E1, float *E2, float *E3);
  int RadiationSource(EnzoVector *extsrc, float *time, float *FS_NGammaDot);
  int GasEnergySource(EnzoVector *extsrc, float *time);
  int ChemistrySource(EnzoVector *extsrc, float *time);
  int AnalyticInitGuess(EnzoVector *u);
  int AnalyticChemistry(EnzoVector *sol, EnzoVector *src, float dt);
  int LinearSolve(EnzoVector *sol, EnzoVector *extsrc, float dt);
  int FSRadiationMask(float *Ef, float *E1, float *E2, float *E3);


 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int BdryType[3][2];

  ///////////////////////////////////////
  // Module-Specific Routines

  // Constructor
  MFSplit();
  
  // Destructor
  ~MFSplit();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem setup
//  int Evolve(HierarchyEntry *ThisGrid, float deltat);
//  int Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat);
  int Evolve(LevelHierarchyEntry *LevelArray[], int level, 
	     HierarchyEntry *Grids[], int NumberOfGrids,
	     TopGridData *MetaData, ExternalBoundary *Exterior, 
#ifdef FAST_SIB
	     SiblingGridList SiblingList[],
#endif
	     float deltat);

  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem debug (for output upon failure)
  int Dump(EnzoVector *ucur);
  
  // Radiation Initializer
  int RadInit(float *Ef, float *E1, float *E2, float *E3, float *E1Units, 
	      float *E2Units, float *E3Units, float *chiint, float *chinuint);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Fill in initial guess for time-evolved solution
  int InitialGuess(EnzoVector *uvec);
  
  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(EnzoVector *unew, int flag);

};


#endif
#endif
