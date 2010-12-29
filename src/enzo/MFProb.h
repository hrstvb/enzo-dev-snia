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
/  Multi-Frequency, Multi-species, Implicit Problem Class
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
#ifndef MF_IMPLICIT_PROBLEM_DEFINED__
#define MF_IMPLICIT_PROBLEM_DEFINED__

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
#include "NonlinearProblemABC.h"
#include "FSProb.h"


class MFProb : public virtual NonlinearProblemABC {

 private:
  
  // flag denoted problem preparedness
  bool prepared;
  
  // overall time spent in solver
  float RTtime;
  float timers[30];
  
  // Free-streaming problem solver, old/new time steps
  FSProb *FSSolve;
  float *Efree;

  // HYPRE Struct-specific data
  Eint32 stSize;                 // stencil size (per frequency)
#ifdef USE_HYPRE
  HYPRE_SStructGrid grid;        // HYPRE grid object for setup
  HYPRE_SStructGraph graph;      // HYPRE graph object for setup
  HYPRE_SStructStencil stencil1; // E1 stencil object
  HYPRE_SStructStencil stencil2; // E2 stencil object
  HYPRE_SStructStencil stencil3; // E3 stencil object
#endif

  // HYPRE Solver-specific data
  Eint32 sol_zeroguess;          // use a zero initial guess
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_relch;              // relative change stopping criteria
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
  HYPRE_SStructMatrix P;         // holds Schur complement matrix
  HYPRE_SStructVector rhsvec;    // holds Schur complement rhs vector
  HYPRE_SStructVector solvec;    // holds Schur complement solution vector
#endif
  Eflt64 *P11tmpvec;             // holds E1 Schur complement matrix entries wrt E1
  Eflt64 *P12tmpvec;             // holds E1 Schur complement matrix entries wrt E2
  Eflt64 *P13tmpvec;             // holds E1 Schur complement matrix entries wrt E3
  Eflt64 *P21tmpvec;             // holds E2 Schur complement matrix entries wrt E1
  Eflt64 *P22tmpvec;             // holds E2 Schur complement matrix entries wrt E2
  Eflt64 *P23tmpvec;             // holds E2 Schur complement matrix entries wrt E3
  Eflt64 *P31tmpvec;             // holds E3 Schur complement matrix entries wrt E1
  Eflt64 *P32tmpvec;             // holds E3 Schur complement matrix entries wrt E2
  Eflt64 *P33tmpvec;             // holds E3 Schur complement matrix entries wrt E3
  Eflt64 *r1tmpvec;              // holds E1 Schur complement rhs entries
  Eflt64 *r2tmpvec;              // holds E2 Schur complement rhs entries
  Eflt64 *r3tmpvec;              // holds E3 Schur complement rhs entries
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // Inexact Newton solver-specific data
  InexactNewtonSolver *INSolve;  // InexactNewton solver
  int approx_jac;                // param to approximate the local jacobian: 
                                 //    0 -> use the analytical jac (not implemented yet)
                                 //    1 -> approximate the jac (default)
  int initial_guess;             // parameter for setting the initial guess:
                                 //    0 -> use previous time step
                                 //    1 -> full fwd Euler (rad) + analytic (other)
                                 //    2 -> analytic (all)
  int semi_implicit;             // use a semi-implicit solver (fixed # of
                                 // Newton iters) as opposed to a fully
                                 // consistent nonlinear solver
                                 //    0 -> use full Newton solver
                                 //  !=0 -> use a fixed # of Newton iters
  int newt_linesearch;           // use a linesearch in the Newton method
                                 //    0 -> no linesearch
                                 //  !=0 -> linesearch (i.e. damped Newton)
  int newt_maxit;                // maximum number of iterations
  int newt_norm;                 // norm for convergence measurement
  float newt_INconst;            // Inexact-Newton constant
  float newt_tol;                // Newton tolerance
  float newt_MinLinesearch;      // minimum allowed line-search length

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of self-gravity problem
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
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  int LimImp;          // implicitness of flux limiter:
                       //    0 -> fully lagged to previous time step
                       //    1 -> fully lagged to previous newton iterate
                       //    2 -> lag only temperature dependence
  int LimType;         // flux limiter formulation:
                       //    0 -> standard Levermore-Pomraning limiter (LP, 1981)
                       //    1 -> rational approx. to LP limiter (LP, 1981)
                       //    2 -> Reynolds approx. to LP limiter
                       //    3 -> no limiter
                       //    4 -> ZEUS limiter (like 1, but no 'albedo')
  EnzoVector *sol;     // solution vector
  EnzoVector *U0;      // old time-level state
  EnzoVector *extsrc;  // temporary vector holding external forcing sources
  EnzoVector *tmp1;    // temporary (see if needed)
  EnzoVector *tmp2;    // temporary (see if needed)
  EnzoVector *tmp3;    // temporary (see if needed)

  int AnalyticChem;    // use analytical reaction solver instead of theta-method
  
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
  float IonizationParms[5];  // parameters for configuring ionization problems

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  float rScale;        // radiation energy density scaling factor
  float eScale;        // specific energy correction scaling factor
  float nScale;        // species density scaling factor
  float EiScale;       // monochromatic->grey radiation scaling (based on spectrum)
  float fsUnits;       // free-streaming radiation unit conversion factor
  float rUnits;        // radiation energy density unit conversion factor
  float rUnits0;       // radiation energy density unit conversion factor (old time)
  float eUnits;        // specific energy correction unit conversion factor
  float nUnits;        // species density unit conversion factor
  float nUnits0;       // species density unit conversion factor (old time)

  float DenUnits;      // density scaling factor
  float DenUnits0;     // density scaling factor (old time)
  float LenUnits;      // length scaling factor
  float LenUnits0;     // length scaling factor (old time)
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor

  // chemistry constants
  float HFrac;         // Fraction of matter composed of Hydrogen

  // storage for integrals over radiation spectrum
  int ESpectrum;       // radiation spectrum type 
                       // (see MFProb_ComputeRadiationIntegrals.src90)
  float hnu0_HI;       // HI ionization threshold (eV)
  float hnu0_HeI;      // HeI ionization threshold (eV)
  float hnu0_HeII;     // HeII ionization threshold (eV)
  float *piHI;         // HI photo-ionization field
  float *piHeI;        // HeI photo-ionization field
  float *piHeII;       // HeII photo-ionization field
  float *GHI;          // HI photo-heating coefficient field
  float *GHeI;         // HeI photo-heating coefficient field
  float *GHeII;        // HeII photo-heating coefficient field

  // linear solver/Jacobian temporary arrays
  EnzoVector *L_e;     // local Jacobian components for energy eqn
  EnzoVector *L_HI;    // local Jacobian components for HI eqn
  EnzoVector *L_HeI;   // local Jacobian components for HeI eqn
  EnzoVector *L_HeII;  // local Jacobian components for HeII eqn
  float *L_E1_E1;      // local Jacobian of E1 eqn wrt E1
  float *L_E2_E2;      // local Jacobian of E2 eqn wrt E2
  float *L_E3_E3;      // local Jacobian of E3 eqn wrt E3
  float *L_E1_HI;      // local Jacobian of E1 eqn wrt HI
  float *L_E2_HI;      // local Jacobian of E2 eqn wrt HI
  float *L_E2_HeI;     // local Jacobian of E2 eqn wrt HeI
  float *L_E3_HI;      // local Jacobian of E3 eqn wrt HI
  float *L_E3_HeI;     // local Jacobian of E3 eqn wrt HeI
  float *L_E3_HeII;    // local Jacobian of E3 eqn wrt HeII

  // access to Enzo data
  float *vx;           // x0-directional velocity
  float *vy;           // x1-directional velocity
  float *vz;           // x2-directional velocity
  float *rho;          // density
  float *eh;           // fluid energy (total or internal)

  // stored arrays for increased efficiency
  float *eCorr;        // gas energy correction
  
  
  // private computation routines
  int LocResid(EnzoVector *locresid, EnzoVector *u);
  int SetupSystem(Eflt64 *m11, Eflt64 *m12, Eflt64 *m13, Eflt64 *m21, 
		  Eflt64 *m22, Eflt64 *m23, Eflt64 *m31, Eflt64 *m32, 
		  Eflt64 *m33, Eflt64 *r1, Eflt64 *r2, Eflt64 *r3, 
		  EnzoVector *u, float *adj11, float *adj12, float *adj13, 
		  float *adj21, float *adj22, float *adj23, float *adj31, 
		  float *adj32, float *adj33);
  int ComputeRadiationIntegrals(EnzoVector *u);
  int EnforceRadiationBounds(float *Ef, float *E1, float *E2, float *E3);
  int RadResid(EnzoVector *radresid, EnzoVector *u);
  int RadJac(EnzoVector *u);
  int BlockSolve(float *Amat, float *xvec, float *bvec, int *N, int *M);
  int Sources(EnzoVector *extsrc, float *time, float *FS_NGammaDot);
  int AnalyticInitGuess(EnzoVector *u);
  int AnalyticResid(EnzoVector *fu, EnzoVector *u);


 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int BdryType[3][2];

  ///////////////////////////////////////
  // Module-Specific Routines

  // Constructor
  MFProb();
  
  // Destructor
  ~MFProb();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem setup
//  int Evolve(HierarchyEntry *ThisGrid, float deltat);
  int Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat);

  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem debug (for output upon failure)
  int Dump(EnzoVector *ucur);
  
  // Radiation Initializer
  int RadInit(float *Ef, float *E1, float *E2, float *E3, float *EiScale);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Enforce boundary conditions onto a vector
  int EnforceBoundary(EnzoVector *vec, int flag);

  // Fill in initial guess for time-evolved solution
  int InitialGuess(EnzoVector *uvec);
  
  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(EnzoVector *unew);

  ///////////////////////////////////////
  // Implicit Solver Interface Routines
  
  // Problem-defining nonlinear residual operations (called repeatedly)
  int nlresid(EnzoVector *fu, EnzoVector *u);
  
  // Problem-specific Linear system setup function, sets up the 
  //   linear Newton system matrix J(u) given an updated state u
  //   (called once per Newton iteration)
  int lsetup(EnzoVector *u);
  
  // Problem-specific Linear solver function 
  //   solves J(u)*s = b to tolerance delta
  //   (called once per Newton iteration)
  int lsolve(EnzoVector *s, EnzoVector *b, EnzoVector *u, float delta);
  
};


#endif
#endif
