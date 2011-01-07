/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit gray flux-limited diffusion solve.
/
/           The variables are stored in the following order: 
/              0 -> radiation energy density
/              1 -> fluid energy correction
/              2:Nspecies+1 -> chemical species (Nspecies may be 0)
/
************************************************************************/

#ifdef TRANSFER
#ifndef AMRFLD_SPLIT_PROBLEM_DEFINED__
#define AMRFLD_SPLIT_PROBLEM_DEFINED__

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
#include "ImplicitProblemABC.h"


class AMRFLDSplit : public virtual ImplicitProblemABC {

 private:
  
  // overall time spent in solver and components
  float RTtime;
  float AMRSolTime;
  
  // AMRsolve-specific data
#ifdef AMR_SOLVE
  AMRsolve_Parameters* amrsolve_params;
#endif

  // solver-specific parameters
  float  sol_tolerance;          // desired solver tolerance
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_type;               // HYPRE solver type
                                 //    0 -> FAC
                                 //    1 -> BiCGStab
                                 //    2 -> BiCGStab-BoomerAMG
                                 //    3 -> GMRES
                                 //    4 -> PFMG (unigrid only)
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging

  // FAC/PFMG-specific solver parameters
  Eint32 sol_rlxtype;            // relaxation type (PFMG only):
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps

  // amrsolve diagnostics
  int totIters;                  // cumulative iterations for solves

  // General problem grid information
  bool OnBdry[3][2];       // denotes if proc owns piece of boundary
  int rank;                // Rank of problem
  int LocDims[3];          // top grid local dims (no ghost or bdry cells)
  float *BdryVals[3][2];   // boundary values for radiation BCs

  // time-stepping related data
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation time step size
  float mindt;         // minimum radiation time step size
  float maxsubcycles;  // max subcycle factor for rad time step within hydro step
  float dtfac;         // desired relative change in radiation per step
  float dtnorm;        // norm choice for computing relative change:
                       //    0 -> max pointwise norm (default)
                       //   >0 -> rms p-norm over entire domain
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dtrad;         // radiation time step size (subcycled)
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  
  // problem defining data
  int Nchem;           // number of chemical species (non-negative integer)
  float NGammaDot;     // ionization strength (photons/sec)
  float EtaRadius;     // ionization source radius
  float EtaCenter[3];  // ionization source location

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  float ErScale;       // radiation energy density scaling factor
  float ErUnits;       // radiation energy density unit conversion factor
  float ErUnits0;      // radiation energy density unit conversion factor
  float NiUnits;       // species density unit conversion factor
  float NiUnits0;      // species density unit conversion factor

  float DenUnits;      // density scaling factor
  float DenUnits0;     // density scaling factor
  float LenUnits;      // length scaling factor
  float LenUnits0;     // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor

  // storage for integrals over radiation spectrum (set during initialization)
  float hnu0_HI;            // HI ionization threshold (eV)
  float hnu0_HeI;           // HeI ionization threshold (eV)
  float hnu0_HeII;          // HeII ionization threshold (eV)
  int ESpectrum;            // integer flag determining spectrum choice
                            //   1 -> 1e5 black body spectrum
                            //   0 -> simple power law spectrum
                            //  -1 -> monochromatic spectrum
  float intSigE;            // int_{nu0}^{inf} sigma_E(nu) d nu
  float intSigESigHI;       // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu) d nu
  float intSigESigHeI;      // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu) d nu
  float intSigESigHeII;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu) d nu
  float intSigESigHInu;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu)/nu d nu
  float intSigESigHeInu;    // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu)/nu d nu
  float intSigESigHeIInu;   // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu)/nu d nu

  // private computation routines
  int EnforceBoundary(LevelHierarchyEntry *LevelArray[]);
  float RadiationSource(LevelHierarchyEntry *LevelArray[], int level, float time);
  float RadiationSpectrum(float nu);
  float CrossSections(float nu, int species);
  int ComputeRadiationIntegrals();
  int FillRates(LevelHierarchyEntry *LevelArray[], int level);
#ifdef AMR_SOLVE
  int RadStep(LevelHierarchyEntry *LevelArray[], int level, 
	      AMRsolve_Hypre_FLD amrfldsolve);
#endif

 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  Eint32 BdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  AMRFLDSplit();
  
  // Destructor
  ~AMRFLDSplit();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Evolver
  int Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat);
  
  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem debug (for output upon failure)
  int Dump(EnzoVector *ucur);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(LevelHierarchyEntry *LevelArray[], int level);

};


#endif
#endif
