#
# FSProb Parameter File:
#   Multi-source free streaming test
#
#  Daniel R. Reynolds, reynolds@smu.edu
#
#############################
#
# Problem-specific initialization parameters
# 
FSProbDensity = 1.67262171e-27       // density (desired nH=1e-3)
FSProbTEnergy = 1.2381613772070e10   // energy [ergs/g] (desired 100 K temp)
FSProbRadiationEnergy = 1.0e-30      // radiation energy (near zero)
FSRadiationOpacity    = 1e-31        // opacity
FSRadiationNGammaDot  = 5.0e48       // avg emissivity source strength 
FSRadiationEtaRadius  = 2.0          // # of sources per proc
#
# General module and solver parameters
# 
FSRadiationScaling      = 1.0e-10
FSRadiationTheta        = 1.0        // implicit Euler
FSRadiationLimiterType  = 4          // Zeus limiter
FSRadiationInitialGuess = 0          // use old time step
FSRadiationTolerance    = 1.0e-5     // linear solver tolerance (absolute)
FSRadiationMaxMGIters   = 50         // allowed Multigrid iters
FSRadiationMGRelaxType  = 2          // MG relaxation type (see HYPRE SysPFMG)
FSRadiationMGPreRelax   = 5          // MG pre-relaxation sweeps
FSRadiationMGPostRelax  = 5          // MG post-relaxation sweeps

#############################
