#
# FSProb Parameter File:
#   Point-source free-streaming radiation test
#
#  Daniel R. Reynolds, reynolds@smu.edu
#
#############################
#
# Problem-specific initialization parameters
# 
FSRadiationBoundaryX0Faces = 2 2     // Neumann
FSRadiationBoundaryX1Faces = 2 2
FSRadiationBoundaryX2Faces = 2 2
FSProbDensity = 1.67262171e-27       // density (desired nH=1e-3)
FSProbTEnergy = 1.2381613772070e10   // energy [ergs/g] (desired 100 K)
FSProbRadiationEnergy = 1.0e-30      // radiation energy (near zero)
FSRadiationOpacity    = 1e-31        // opacity
FSRadiationNGammaDot  = 5.0e48       // emissivity source strength
FSRadiationEtaRadius  = 1.0e0        // radius for emissivity source (in cells)
#
# General module and solver parameters
# 
FSRadiationMaxDt        = 0.01       // 1/100 of light-crossing time
FSRadiationScaling      = 1.0e-14    // unit scaling
FSRadiationTheta        = 1.0        // implicit Euler
FSRadiationLimiterType  = 4          // Zeus limiter
FSRadiationInitialGuess = 0          // use old time step
FSRadiationTolerance    = 1.0e-5     // linear solver tolerance (absolute)
FSRadiationMaxMGIters   = 50         // allowed Multigrid iters
FSRadiationMGRelaxType  = 2          // MG relaxation type (see HYPRE PFMG)
FSRadiationMGPreRelax   = 5          // MG pre-relaxation sweeps
FSRadiationMGPostRelax  = 5          // MG post-relaxation sweeps

#############################
