
# Automatically-generated global data handler
# DO NOT MODIFY DIRECTLY

cdef extern from "units.h":
    pass

cdef extern from "flowdefs.h":
    pass

cdef extern from "CosmologyParameters.h":
    pass

cdef extern from "communication.h":
    pass

cdef extern from "StarParticleData.h":
    pass

cdef extern from "global_data.h":
    cdef extern int LoadBalancing
    cdef extern int LoadBalancingCycleSkip
    cdef extern int CoresPerNode
    cdef extern int PreviousMaxTask
    cdef extern int FileDirectedOutput
    cdef extern int debug
    cdef extern int debug1
    cdef extern int debug2
    cdef extern int extract
    cdef extern int ProblemType
    cdef extern double huge_number
    cdef extern double tiny_number
    cdef extern double Gamma
    cdef extern int PressureFree
    cdef extern int RefineBy
    cdef extern int MaximumRefinementLevel
    cdef extern int MaximumGravityRefinementLevel
    cdef extern int MaximumParticleRefinementLevel
    cdef extern int CellFlaggingMethod[MAX_FLAGGING_METHODS]
    cdef extern double MustRefineRegionLeftEdge[MAX_DIMENSION]
    cdef extern double MustRefineRegionRightEdge[MAX_DIMENSION]
    cdef extern int MustRefineRegionMinRefinementLevel
    cdef extern int MetallicityRefinementMinLevel
    cdef extern double MetallicityRefinementMinMetallicity
    cdef extern double TimestepSafetyVelocity
    cdef extern int FluxCorrection
    cdef extern int ConservativeInterpolation
    cdef extern double MinimumEfficiency
    cdef extern int MinimumSubgridEdge
    cdef extern int MaximumSubgridSize
    cdef extern int NumberOfBufferZones
    cdef extern double DomainLeftEdge[MAX_DIMENSION]
    cdef extern double DomainRightEdge[MAX_DIMENSION]
    cdef extern double GridVelocity[MAX_DIMENSION]
    cdef extern double RefineRegionLeftEdge[MAX_DIMENSION]
    cdef extern double RefineRegionRightEdge[MAX_DIMENSION]
    cdef extern int RefineRegionAutoAdjust
    cdef extern int UniformGravity
    cdef extern int UniformGravityDirection
    cdef extern double UniformGravityConstant
    cdef extern int PointSourceGravity
    cdef extern double PointSourceGravityPosition[MAX_DIMENSION]
    cdef extern double PointSourceGravityConstant
    cdef extern double PointSourceGravityCoreRadius
    cdef extern int SelfGravity
    cdef extern int CopyGravPotential
    cdef extern int PotentialIterations
    cdef extern int BaryonSelfGravityApproximation
    cdef extern double GravitationalConstant
    cdef extern double S2ParticleSize
    cdef extern double GravityResolution
    cdef extern int ComputePotential
    cdef extern int WritePotential
    cdef extern int GreensFunctionMaxNumber
    cdef extern int GreensFunctionMaxSize
    cdef extern int DualEnergyFormalism
    cdef extern double DualEnergyFormalismEta1
    cdef extern double DualEnergyFormalismEta2
    cdef extern double ParticleCourantSafetyNumber
    cdef extern double RootGridCourantSafetyNumber
    cdef extern int RadiativeCooling
    cdef extern int GadgetEquilibriumCooling
    cdef extern int RandomForcing
    cdef extern double RandomForcingEdot
    cdef extern double RandomForcingMachNumber
    cdef extern int MultiSpecies
    cdef extern int GloverChemistryModel
    cdef extern int GloverRadiationBackground
    cdef extern int GloverOpticalDepth
    cdef extern int MultiMetals
    cdef extern int RadiationFieldType
    cdef extern int AdjustUVBackground
    cdef extern double SetUVBAmplitude
    cdef extern double SetHeIIHeatingScale
    cdef extern int RadiationFieldLevelRecompute
    cdef extern int RadiationXRaySecondaryIon
    cdef extern int OutputCoolingTime
    cdef extern int OutputTemperature
    cdef extern int OutputSmoothedDarkMatter
    cdef extern int SmoothedDarkMatterNeighbors
    cdef extern double ZEUSLinearArtificialViscosity
    cdef extern double ZEUSQuadraticArtificialViscosity
    cdef extern int UseMinimumPressureSupport
    cdef extern double MinimumPressureSupportParameter
    cdef extern double StaticRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION]
    cdef extern double StaticRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION]
    cdef extern int StaticRefineRegionLevel[MAX_STATIC_REGIONS]
    cdef extern int MyProcessorNumber
    cdef extern int NumberOfProcessors
    cdef extern double CommunicationTime
    cdef extern int ParallelRootGridIO
    cdef extern int ParallelParticleIO
    cdef extern int Unigrid
    cdef extern int CubeDumpEnabled
    cdef extern int PartitionNestedGrids
    cdef extern int ExtractFieldsOnly
    cdef extern int First_Pass
    cdef extern int UnigridTranspose
    cdef extern int PythonSubcycleSkip
    cdef extern int InlineHaloFinder
    cdef extern int HaloFinderSubfind
    cdef extern int HaloFinderOutputParticleList
    cdef extern int HaloFinderMinimumSize
    cdef extern int HaloFinderCycleSkip
    cdef extern double HaloFinderLinkingLength
    cdef extern double HaloFinderTimestep
    cdef extern double HaloFinderLastTime
    cdef extern double MinimumSlopeForRefinement[MAX_FLAGGING_METHODS]
    cdef extern int SlopeFlaggingFields[MAX_FLAGGING_METHODS]
    cdef extern double MinimumOverDensityForRefinement[MAX_FLAGGING_METHODS]
    cdef extern double MinimumMassForRefinement[MAX_FLAGGING_METHODS]
    cdef extern double MinimumMassForRefinementLevelExponent[MAX_FLAGGING_METHODS]
    cdef extern double DepositPositionsParticleSmoothRadius
    cdef extern double MinimumPressureJumpForRefinement
    cdef extern double MinimumEnergyRatioForRefinement
    cdef extern double RefineByJeansLengthSafetyFactor
    cdef extern int MustRefineParticlesRefineToLevel
    cdef extern double MinimumShearForRefinement
    cdef extern double RefineByResistiveLengthSafetyFactor
    cdef extern int NohProblemFullBox
    cdef extern int ComovingCoordinates
    cdef extern int StarParticleCreation
    cdef extern int StarParticleFeedback
    cdef extern int NumberOfParticleAttributes
    cdef extern int AddParticleAttributes
    cdef extern int TimeActionType[MAX_TIME_ACTIONS]
    cdef extern double TimeActionTime[MAX_TIME_ACTIONS]
    cdef extern double TimeActionRedshift[MAX_TIME_ACTIONS]
    cdef extern double TimeActionParameter[MAX_TIME_ACTIONS]
    cdef extern int TracerParticleOn
    cdef extern double TracerParticleCreationSpacing
    cdef extern double TracerParticleCreationLeftEdge[MAX_DIMENSION]
    cdef extern double TracerParticleCreationRightEdge[MAX_DIMENSION]
    cdef extern int ParticleTypeInFile
    cdef extern int OutputParticleTypeGrouping
    cdef extern int ExternalBoundaryIO
    cdef extern int ExternalBoundaryTypeIO
    cdef extern int ExternalBoundaryValueIO
    cdef extern int ExternalBoundaryField
    cdef extern int SimpleConstantBoundary
    cdef extern int TaskMap[MAX_NUMBER_OF_TASKS]
    cdef extern double NodeMem[MAX_NUMBER_OF_NODES]
    cdef extern int NodeMap[MAX_NUMBER_OF_NODES]
    cdef extern int LoadGridDataAtStart
    cdef extern double GlobalCommunication
    cdef extern double RecvComm
    cdef extern double WaitComm
    cdef extern double timer[MAX_COUNTERS]
    cdef extern int counter[MAX_COUNTERS]
    cdef extern int traceMEM
    cdef extern double starttime
    cdef extern double endtime
    cdef extern double Start_Wall_Time
    cdef extern double End_Wall_Time
    cdef extern double WallTime
    cdef extern int flagging_count
    cdef extern int in_count
    cdef extern int out_count
    cdef extern int moving_count
    cdef extern double flagging_pct
    cdef extern double moving_pct
    cdef extern int traceMPI
    cdef extern int MovieDataField[MAX_MOVIE_FIELDS]
    cdef extern int MovieSkipTimestep
    cdef extern int Movie3DVolumes
    cdef extern int MovieVertexCentered
    cdef extern int NewMovieDumpNumber
    cdef extern int NewMovieParticleOn
    cdef extern int UseHydro
    cdef extern int Coordinate
    cdef extern int NSpecies
    cdef extern int NColor
    cdef extern double Theta_Limiter
    cdef extern int RKOrder
    cdef extern int UsePhysicalUnit
    cdef extern int iden
    cdef extern int ietot
    cdef extern int ivx
    cdef extern int ivy
    cdef extern int ivz
    cdef extern int iBx
    cdef extern int iBy
    cdef extern int iBz
    cdef extern int iPhi
    cdef extern int ieint
    cdef extern int iD
    cdef extern int iEtot
    cdef extern int iS1
    cdef extern int iS2
    cdef extern int iS3
    cdef extern int iEint
    cdef extern double SmallRho
    cdef extern double SmallP
    cdef extern double SmallEint
    cdef extern double SmallT
    cdef extern double MaximumAlvenSpeed
    cdef extern int NEQ_HYDRO
    cdef extern int NEQ_MHD
    cdef extern int ReconstructionMethod
    cdef extern int RiemannSolver
    cdef extern int EOSType
    cdef extern double EOSSoundSpeed
    cdef extern double EOSCriticalDensity
    cdef extern double EOSGamma
    cdef extern double C_h
    cdef extern double C_p
    cdef extern int UseConstantAcceleration
    cdef extern double ConstantAcceleration[MAX_DIMENSION]
    cdef extern double Mu
    cdef extern int ExternalGravity
    cdef extern int StringKick
    cdef extern int UseFloor
    cdef extern int UseViscosity
    cdef extern int UseAmbipolarDiffusion
    cdef extern int UseResistivity
    cdef extern int UseH2OnDust
    cdef extern double PhotoelectricHeating
    cdef extern double CoolingCutOffDensity1
    cdef extern double CoolingCutOffDensity2
    cdef extern double CoolingPowerCutOffDensity1
    cdef extern double CoolingPowerCutOffDensity2
    cdef extern double CoolingCutOffTemperature
    cdef extern int CoolingModel
    cdef extern double HaloMass
    cdef extern double HaloConcentration
    cdef extern double HaloRedshift
    cdef extern double HaloCentralDensity
    cdef extern double HaloVirialRadius
    cdef extern double ExternalGravityDensity
    cdef extern double ExternalGravityRadius
    cdef extern int UseDivergenceCleaning
    cdef extern int DivergenceCleaningBoundaryBuffer
    cdef extern double DivergenceCleaningThreshold
    cdef extern double PoissonApproximationThreshold
    cdef extern int ShiningParticleID
    cdef extern double SinkMergeDistance
    cdef extern double SinkMergeMass
    cdef extern double TotalSinkMass
    cdef extern int StellarWindFeedback
    cdef extern double StellarWindTurnOnMass
    cdef extern int NBodyDirectSummation
    cdef extern int UseDrivingField
    cdef extern double DrivingEfficiency
    cdef extern int UseCUDA
    cdef extern int ran1_init
    cdef extern int MemoryLimit
    cdef extern int MetalCooling
    cdef extern int CIECooling
    cdef extern int H2OpticalDepthApproximation
    cdef extern int RadiativeTransfer
    cdef extern double PhotonTime
    cdef extern double dtPhoton
    cdef extern double RadiativeTransferSourceRadius
    cdef extern double RadiativeTransferPropagationSpeedFraction
    cdef extern double RadiativeTransferPropagationDistance
    cdef extern double RadiativeTransferSplitPhotonRadius
    cdef extern double RadiativeTransferRaysPerCell
    cdef extern int RadiativeTransferInitialHEALPixLevel
    cdef extern double RadiativeTransferPhotonEscapeRadius
    cdef extern int RadiativeTransferInterpolateField
    cdef extern double RadiativeTransferTimestepVelocityLimit
    cdef extern int RadiativeTransferSourceClustering
    cdef extern double RadiativeTransferPhotonMergeRadius
    cdef extern int RadiationPressure
    cdef extern int RadiativeTransferOpticallyThinH2
    cdef extern int RadiativeTransferPeriodicBoundary
    cdef extern double EscapedPhotonCount[4]
    cdef extern double TotalEscapedPhotonCount[4]
    cdef extern int FieldsToInterpolate[MAX_NUMBER_OF_BARYON_FIELDS]
    cdef extern int RadiativeTransferCoupledRateSolver
    cdef extern double AngularVelocity
    cdef extern double VelocityGradient
    cdef extern int ShearingBoundaryDirection
    cdef extern int ShearingVelocityDirection
    cdef extern int ShearingOtherDirection
    cdef extern int useMHD
    cdef extern double TopGridDx[MAX_DIMENSION]
    cdef extern int ShearingBoxProblemType
    cdef extern double IsothermalSoundSpeed
    cdef extern int RefineByJeansLengthUnits
    cdef extern double GlobalLengthUnits
    cdef extern double GlobalMassUnits
    cdef extern double GlobalDensityUnits
    cdef extern double GlobalTimeUnits
    cdef extern int flow_trace_on
    cdef extern int flow_trace_level
    cdef extern int CommunicationDirection
    cdef extern double HubbleConstantNow
    cdef extern double OmegaMatterNow
    cdef extern double OmegaLambdaNow
    cdef extern double ComovingBoxSize
    cdef extern double MaxExpansionRate
    cdef extern double InitialTimeInCodeUnits
    cdef extern double InitialRedshift
    cdef extern double CosmologyOutputRedshift[MAX_NUMBER_OF_OUTPUT_REDSHIFTS]
    cdef extern double CosmologyOutputRedshiftTime[MAX_NUMBER_OF_OUTPUT_REDSHIFTS]
    cdef extern int NumberOfStarParticles
    cdef extern int TotalNumberOfStars
    cdef extern double StarMakerOverDensityThreshold
    cdef extern double StarMakerMassEfficiency
    cdef extern double StarMakerMinimumMass
    cdef extern double StarMakerMinimumDynamicalTime
    cdef extern double StarMassEjectionFraction
    cdef extern double StarMetalYield
    cdef extern double StarEnergyToThermalFeedback
    cdef extern double StarEnergyToStellarUV
    cdef extern double StarEnergyToQuasarUV
    cdef extern double PopIIIStarMass
    cdef extern int PopIIIBlackHoles
    cdef extern double PopIIIBHLuminosityEfficiency
    cdef extern double PopIIIOverDensityThreshold
    cdef extern double PopIIIH2CriticalFraction
    cdef extern double PopIIIMetalCriticalFraction
    cdef extern double PopIIISupernovaRadius
    cdef extern int PopIIISupernovaUseColour
    cdef extern int StarClusterUseMetalField
    cdef extern double StarClusterMinDynamicalTime
    cdef extern double StarClusterIonizingLuminosity
    cdef extern double StarClusterSNEnergy
    cdef extern double StarClusterSNRadius
    cdef extern double StarClusterFormEfficiency
    cdef extern double StarClusterMinimumMass
    cdef extern double StarClusterCombineRadius
    cdef extern double StarClusterRegionLeftEdge[MAX_DIMENSION]
    cdef extern double StarClusterRegionRightEdge[MAX_DIMENSION]
    cdef extern double MBHMinDynamicalTime
    cdef extern double MBHMinimumMass
    cdef extern int MBHFeedbackThermal
    cdef extern double MBHFeedbackRadius
    cdef extern double MBHFeedbackRadiativeEfficiency
    cdef extern double MBHFeedbackThermalCoupling
    cdef extern double MBHCombineRadius
    cdef extern double MBHIonizingLuminosity
    cdef extern double minStarLifetime
    cdef extern double LastSupernovaTime
    cdef extern int CommunicationReceiveCurrentDependsOn
    cdef extern int CommunicationReceiveIndex
    cdef extern int CommunicationReceiveCallType[MAX_RECEIVE_BUFFERS]
    cdef extern int CommunicationReceiveDependsOn[MAX_RECEIVE_BUFFERS]
    cdef extern double CommunicationReceiveArgument[MAX_DIMENSION][MAX_RECEIVE_BUFFERS]
    cdef extern int CommunicationReceiveArgumentInt[MAX_DIMENSION][MAX_RECEIVE_BUFFERS]

    cdef extern char *DataLabel[MAX_NUMBER_OF_BARYON_FIELDS]
    cdef extern char *DataUnits[MAX_NUMBER_OF_BARYON_FIELDS]

cdef class LabelHandler:
    cdef char **data
    cdef int number

    def __cinit__(self, int number):
        self.number = number

    def __getitem__(self, int index):
        if index > self.number: raise IndexError
        elif index < 0: raise IndexError
        return <char *> self.data[index]

cdef class _global_data:
    cdef public LabelHandler DataLabel
    cdef public LabelHandler DataUnits

    def __cinit__(self):
        global DataLabel
        global DataUnits
        self.DataLabel = LabelHandler(MAX_NUMBER_OF_BARYON_FIELDS)
        self.DataLabel.data = DataLabel
        self.DataUnits = LabelHandler(MAX_NUMBER_OF_BARYON_FIELDS)
        self.DataUnits.data = DataUnits

    property LoadBalancing:
        def __get__(self):
            global LoadBalancing
            return LoadBalancing
        def __set__(self, val):
            global LoadBalancing
            LoadBalancing = val

    property LoadBalancingCycleSkip:
        def __get__(self):
            global LoadBalancingCycleSkip
            return LoadBalancingCycleSkip
        def __set__(self, val):
            global LoadBalancingCycleSkip
            LoadBalancingCycleSkip = val

    property CoresPerNode:
        def __get__(self):
            global CoresPerNode
            return CoresPerNode
        def __set__(self, val):
            global CoresPerNode
            CoresPerNode = val

    property PreviousMaxTask:
        def __get__(self):
            global PreviousMaxTask
            return PreviousMaxTask
        def __set__(self, val):
            global PreviousMaxTask
            PreviousMaxTask = val

    property FileDirectedOutput:
        def __get__(self):
            global FileDirectedOutput
            return FileDirectedOutput
        def __set__(self, val):
            global FileDirectedOutput
            FileDirectedOutput = val

    property debug:
        def __get__(self):
            global debug
            return debug
        def __set__(self, val):
            global debug
            debug = val

    property debug1:
        def __get__(self):
            global debug1
            return debug1
        def __set__(self, val):
            global debug1
            debug1 = val

    property debug2:
        def __get__(self):
            global debug2
            return debug2
        def __set__(self, val):
            global debug2
            debug2 = val

    property extract:
        def __get__(self):
            global extract
            return extract
        def __set__(self, val):
            global extract
            extract = val

    property ProblemType:
        def __get__(self):
            global ProblemType
            return ProblemType
        def __set__(self, val):
            global ProblemType
            ProblemType = val

    property huge_number:
        def __get__(self):
            global huge_number
            return huge_number
        def __set__(self, val):
            global huge_number
            huge_number = val

    property tiny_number:
        def __get__(self):
            global tiny_number
            return tiny_number
        def __set__(self, val):
            global tiny_number
            tiny_number = val

    property Gamma:
        def __get__(self):
            global Gamma
            return Gamma
        def __set__(self, val):
            global Gamma
            Gamma = val

    property PressureFree:
        def __get__(self):
            global PressureFree
            return PressureFree
        def __set__(self, val):
            global PressureFree
            PressureFree = val

    property RefineBy:
        def __get__(self):
            global RefineBy
            return RefineBy
        def __set__(self, val):
            global RefineBy
            RefineBy = val

    property MaximumRefinementLevel:
        def __get__(self):
            global MaximumRefinementLevel
            return MaximumRefinementLevel
        def __set__(self, val):
            global MaximumRefinementLevel
            MaximumRefinementLevel = val

    property MaximumGravityRefinementLevel:
        def __get__(self):
            global MaximumGravityRefinementLevel
            return MaximumGravityRefinementLevel
        def __set__(self, val):
            global MaximumGravityRefinementLevel
            MaximumGravityRefinementLevel = val

    property MaximumParticleRefinementLevel:
        def __get__(self):
            global MaximumParticleRefinementLevel
            return MaximumParticleRefinementLevel
        def __set__(self, val):
            global MaximumParticleRefinementLevel
            MaximumParticleRefinementLevel = val

    property CellFlaggingMethod:
        def __get__(self):
            print "Returning a copy of CellFlaggingMethod"
            global CellFlaggingMethod
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(CellFlaggingMethod[i])
            return retval
        def __set__(self, val):
            cdef int i
            global CellFlaggingMethod
            for i in range(MAX_FLAGGING_METHODS):
                CellFlaggingMethod[i] = val[i]

    property MustRefineRegionLeftEdge:
        def __get__(self):
            print "Returning a copy of MustRefineRegionLeftEdge"
            global MustRefineRegionLeftEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(MustRefineRegionLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MustRefineRegionLeftEdge
            for i in range(MAX_DIMENSION):
                MustRefineRegionLeftEdge[i] = val[i]

    property MustRefineRegionRightEdge:
        def __get__(self):
            print "Returning a copy of MustRefineRegionRightEdge"
            global MustRefineRegionRightEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(MustRefineRegionRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MustRefineRegionRightEdge
            for i in range(MAX_DIMENSION):
                MustRefineRegionRightEdge[i] = val[i]

    property MustRefineRegionMinRefinementLevel:
        def __get__(self):
            global MustRefineRegionMinRefinementLevel
            return MustRefineRegionMinRefinementLevel
        def __set__(self, val):
            global MustRefineRegionMinRefinementLevel
            MustRefineRegionMinRefinementLevel = val

    property MetallicityRefinementMinLevel:
        def __get__(self):
            global MetallicityRefinementMinLevel
            return MetallicityRefinementMinLevel
        def __set__(self, val):
            global MetallicityRefinementMinLevel
            MetallicityRefinementMinLevel = val

    property MetallicityRefinementMinMetallicity:
        def __get__(self):
            global MetallicityRefinementMinMetallicity
            return MetallicityRefinementMinMetallicity
        def __set__(self, val):
            global MetallicityRefinementMinMetallicity
            MetallicityRefinementMinMetallicity = val

    property TimestepSafetyVelocity:
        def __get__(self):
            global TimestepSafetyVelocity
            return TimestepSafetyVelocity
        def __set__(self, val):
            global TimestepSafetyVelocity
            TimestepSafetyVelocity = val

    property FluxCorrection:
        def __get__(self):
            global FluxCorrection
            return FluxCorrection
        def __set__(self, val):
            global FluxCorrection
            FluxCorrection = val

    property ConservativeInterpolation:
        def __get__(self):
            global ConservativeInterpolation
            return ConservativeInterpolation
        def __set__(self, val):
            global ConservativeInterpolation
            ConservativeInterpolation = val

    property MinimumEfficiency:
        def __get__(self):
            global MinimumEfficiency
            return MinimumEfficiency
        def __set__(self, val):
            global MinimumEfficiency
            MinimumEfficiency = val

    property MinimumSubgridEdge:
        def __get__(self):
            global MinimumSubgridEdge
            return MinimumSubgridEdge
        def __set__(self, val):
            global MinimumSubgridEdge
            MinimumSubgridEdge = val

    property MaximumSubgridSize:
        def __get__(self):
            global MaximumSubgridSize
            return MaximumSubgridSize
        def __set__(self, val):
            global MaximumSubgridSize
            MaximumSubgridSize = val

    property NumberOfBufferZones:
        def __get__(self):
            global NumberOfBufferZones
            return NumberOfBufferZones
        def __set__(self, val):
            global NumberOfBufferZones
            NumberOfBufferZones = val

    property DomainLeftEdge:
        def __get__(self):
            print "Returning a copy of DomainLeftEdge"
            global DomainLeftEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(DomainLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global DomainLeftEdge
            for i in range(MAX_DIMENSION):
                DomainLeftEdge[i] = val[i]

    property DomainRightEdge:
        def __get__(self):
            print "Returning a copy of DomainRightEdge"
            global DomainRightEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(DomainRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global DomainRightEdge
            for i in range(MAX_DIMENSION):
                DomainRightEdge[i] = val[i]

    property GridVelocity:
        def __get__(self):
            print "Returning a copy of GridVelocity"
            global GridVelocity
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(GridVelocity[i])
            return retval
        def __set__(self, val):
            cdef int i
            global GridVelocity
            for i in range(MAX_DIMENSION):
                GridVelocity[i] = val[i]

    property RefineRegionLeftEdge:
        def __get__(self):
            print "Returning a copy of RefineRegionLeftEdge"
            global RefineRegionLeftEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(RefineRegionLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global RefineRegionLeftEdge
            for i in range(MAX_DIMENSION):
                RefineRegionLeftEdge[i] = val[i]

    property RefineRegionRightEdge:
        def __get__(self):
            print "Returning a copy of RefineRegionRightEdge"
            global RefineRegionRightEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(RefineRegionRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global RefineRegionRightEdge
            for i in range(MAX_DIMENSION):
                RefineRegionRightEdge[i] = val[i]

    property RefineRegionAutoAdjust:
        def __get__(self):
            global RefineRegionAutoAdjust
            return RefineRegionAutoAdjust
        def __set__(self, val):
            global RefineRegionAutoAdjust
            RefineRegionAutoAdjust = val

    property UniformGravity:
        def __get__(self):
            global UniformGravity
            return UniformGravity
        def __set__(self, val):
            global UniformGravity
            UniformGravity = val

    property UniformGravityDirection:
        def __get__(self):
            global UniformGravityDirection
            return UniformGravityDirection
        def __set__(self, val):
            global UniformGravityDirection
            UniformGravityDirection = val

    property UniformGravityConstant:
        def __get__(self):
            global UniformGravityConstant
            return UniformGravityConstant
        def __set__(self, val):
            global UniformGravityConstant
            UniformGravityConstant = val

    property PointSourceGravity:
        def __get__(self):
            global PointSourceGravity
            return PointSourceGravity
        def __set__(self, val):
            global PointSourceGravity
            PointSourceGravity = val

    property PointSourceGravityPosition:
        def __get__(self):
            print "Returning a copy of PointSourceGravityPosition"
            global PointSourceGravityPosition
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(PointSourceGravityPosition[i])
            return retval
        def __set__(self, val):
            cdef int i
            global PointSourceGravityPosition
            for i in range(MAX_DIMENSION):
                PointSourceGravityPosition[i] = val[i]

    property PointSourceGravityConstant:
        def __get__(self):
            global PointSourceGravityConstant
            return PointSourceGravityConstant
        def __set__(self, val):
            global PointSourceGravityConstant
            PointSourceGravityConstant = val

    property PointSourceGravityCoreRadius:
        def __get__(self):
            global PointSourceGravityCoreRadius
            return PointSourceGravityCoreRadius
        def __set__(self, val):
            global PointSourceGravityCoreRadius
            PointSourceGravityCoreRadius = val

    property SelfGravity:
        def __get__(self):
            global SelfGravity
            return SelfGravity
        def __set__(self, val):
            global SelfGravity
            SelfGravity = val

    property CopyGravPotential:
        def __get__(self):
            global CopyGravPotential
            return CopyGravPotential
        def __set__(self, val):
            global CopyGravPotential
            CopyGravPotential = val

    property PotentialIterations:
        def __get__(self):
            global PotentialIterations
            return PotentialIterations
        def __set__(self, val):
            global PotentialIterations
            PotentialIterations = val

    property BaryonSelfGravityApproximation:
        def __get__(self):
            global BaryonSelfGravityApproximation
            return BaryonSelfGravityApproximation
        def __set__(self, val):
            global BaryonSelfGravityApproximation
            BaryonSelfGravityApproximation = val

    property GravitationalConstant:
        def __get__(self):
            global GravitationalConstant
            return GravitationalConstant
        def __set__(self, val):
            global GravitationalConstant
            GravitationalConstant = val

    property S2ParticleSize:
        def __get__(self):
            global S2ParticleSize
            return S2ParticleSize
        def __set__(self, val):
            global S2ParticleSize
            S2ParticleSize = val

    property GravityResolution:
        def __get__(self):
            global GravityResolution
            return GravityResolution
        def __set__(self, val):
            global GravityResolution
            GravityResolution = val

    property ComputePotential:
        def __get__(self):
            global ComputePotential
            return ComputePotential
        def __set__(self, val):
            global ComputePotential
            ComputePotential = val

    property WritePotential:
        def __get__(self):
            global WritePotential
            return WritePotential
        def __set__(self, val):
            global WritePotential
            WritePotential = val

    property GreensFunctionMaxNumber:
        def __get__(self):
            global GreensFunctionMaxNumber
            return GreensFunctionMaxNumber
        def __set__(self, val):
            global GreensFunctionMaxNumber
            GreensFunctionMaxNumber = val

    property GreensFunctionMaxSize:
        def __get__(self):
            global GreensFunctionMaxSize
            return GreensFunctionMaxSize
        def __set__(self, val):
            global GreensFunctionMaxSize
            GreensFunctionMaxSize = val

    property DualEnergyFormalism:
        def __get__(self):
            global DualEnergyFormalism
            return DualEnergyFormalism
        def __set__(self, val):
            global DualEnergyFormalism
            DualEnergyFormalism = val

    property DualEnergyFormalismEta1:
        def __get__(self):
            global DualEnergyFormalismEta1
            return DualEnergyFormalismEta1
        def __set__(self, val):
            global DualEnergyFormalismEta1
            DualEnergyFormalismEta1 = val

    property DualEnergyFormalismEta2:
        def __get__(self):
            global DualEnergyFormalismEta2
            return DualEnergyFormalismEta2
        def __set__(self, val):
            global DualEnergyFormalismEta2
            DualEnergyFormalismEta2 = val

    property ParticleCourantSafetyNumber:
        def __get__(self):
            global ParticleCourantSafetyNumber
            return ParticleCourantSafetyNumber
        def __set__(self, val):
            global ParticleCourantSafetyNumber
            ParticleCourantSafetyNumber = val

    property RootGridCourantSafetyNumber:
        def __get__(self):
            global RootGridCourantSafetyNumber
            return RootGridCourantSafetyNumber
        def __set__(self, val):
            global RootGridCourantSafetyNumber
            RootGridCourantSafetyNumber = val

    property RadiativeCooling:
        def __get__(self):
            global RadiativeCooling
            return RadiativeCooling
        def __set__(self, val):
            global RadiativeCooling
            RadiativeCooling = val

    property GadgetEquilibriumCooling:
        def __get__(self):
            global GadgetEquilibriumCooling
            return GadgetEquilibriumCooling
        def __set__(self, val):
            global GadgetEquilibriumCooling
            GadgetEquilibriumCooling = val

    property RandomForcing:
        def __get__(self):
            global RandomForcing
            return RandomForcing
        def __set__(self, val):
            global RandomForcing
            RandomForcing = val

    property RandomForcingEdot:
        def __get__(self):
            global RandomForcingEdot
            return RandomForcingEdot
        def __set__(self, val):
            global RandomForcingEdot
            RandomForcingEdot = val

    property RandomForcingMachNumber:
        def __get__(self):
            global RandomForcingMachNumber
            return RandomForcingMachNumber
        def __set__(self, val):
            global RandomForcingMachNumber
            RandomForcingMachNumber = val

    property MultiSpecies:
        def __get__(self):
            global MultiSpecies
            return MultiSpecies
        def __set__(self, val):
            global MultiSpecies
            MultiSpecies = val

    property GloverChemistryModel:
        def __get__(self):
            global GloverChemistryModel
            return GloverChemistryModel
        def __set__(self, val):
            global GloverChemistryModel
            GloverChemistryModel = val

    property GloverRadiationBackground:
        def __get__(self):
            global GloverRadiationBackground
            return GloverRadiationBackground
        def __set__(self, val):
            global GloverRadiationBackground
            GloverRadiationBackground = val

    property GloverOpticalDepth:
        def __get__(self):
            global GloverOpticalDepth
            return GloverOpticalDepth
        def __set__(self, val):
            global GloverOpticalDepth
            GloverOpticalDepth = val

    property MultiMetals:
        def __get__(self):
            global MultiMetals
            return MultiMetals
        def __set__(self, val):
            global MultiMetals
            MultiMetals = val

    property RadiationFieldType:
        def __get__(self):
            global RadiationFieldType
            return RadiationFieldType
        def __set__(self, val):
            global RadiationFieldType
            RadiationFieldType = val

    property AdjustUVBackground:
        def __get__(self):
            global AdjustUVBackground
            return AdjustUVBackground
        def __set__(self, val):
            global AdjustUVBackground
            AdjustUVBackground = val

    property SetUVBAmplitude:
        def __get__(self):
            global SetUVBAmplitude
            return SetUVBAmplitude
        def __set__(self, val):
            global SetUVBAmplitude
            SetUVBAmplitude = val

    property SetHeIIHeatingScale:
        def __get__(self):
            global SetHeIIHeatingScale
            return SetHeIIHeatingScale
        def __set__(self, val):
            global SetHeIIHeatingScale
            SetHeIIHeatingScale = val

    property RadiationFieldLevelRecompute:
        def __get__(self):
            global RadiationFieldLevelRecompute
            return RadiationFieldLevelRecompute
        def __set__(self, val):
            global RadiationFieldLevelRecompute
            RadiationFieldLevelRecompute = val

    property RadiationXRaySecondaryIon:
        def __get__(self):
            global RadiationXRaySecondaryIon
            return RadiationXRaySecondaryIon
        def __set__(self, val):
            global RadiationXRaySecondaryIon
            RadiationXRaySecondaryIon = val

    property OutputCoolingTime:
        def __get__(self):
            global OutputCoolingTime
            return OutputCoolingTime
        def __set__(self, val):
            global OutputCoolingTime
            OutputCoolingTime = val

    property OutputTemperature:
        def __get__(self):
            global OutputTemperature
            return OutputTemperature
        def __set__(self, val):
            global OutputTemperature
            OutputTemperature = val

    property OutputSmoothedDarkMatter:
        def __get__(self):
            global OutputSmoothedDarkMatter
            return OutputSmoothedDarkMatter
        def __set__(self, val):
            global OutputSmoothedDarkMatter
            OutputSmoothedDarkMatter = val

    property SmoothedDarkMatterNeighbors:
        def __get__(self):
            global SmoothedDarkMatterNeighbors
            return SmoothedDarkMatterNeighbors
        def __set__(self, val):
            global SmoothedDarkMatterNeighbors
            SmoothedDarkMatterNeighbors = val

    property ZEUSLinearArtificialViscosity:
        def __get__(self):
            global ZEUSLinearArtificialViscosity
            return ZEUSLinearArtificialViscosity
        def __set__(self, val):
            global ZEUSLinearArtificialViscosity
            ZEUSLinearArtificialViscosity = val

    property ZEUSQuadraticArtificialViscosity:
        def __get__(self):
            global ZEUSQuadraticArtificialViscosity
            return ZEUSQuadraticArtificialViscosity
        def __set__(self, val):
            global ZEUSQuadraticArtificialViscosity
            ZEUSQuadraticArtificialViscosity = val

    property UseMinimumPressureSupport:
        def __get__(self):
            global UseMinimumPressureSupport
            return UseMinimumPressureSupport
        def __set__(self, val):
            global UseMinimumPressureSupport
            UseMinimumPressureSupport = val

    property MinimumPressureSupportParameter:
        def __get__(self):
            global MinimumPressureSupportParameter
            return MinimumPressureSupportParameter
        def __set__(self, val):
            global MinimumPressureSupportParameter
            MinimumPressureSupportParameter = val

    property StaticRefineRegionLeftEdge:
        def __get__(self):
            print "Returning a copy of StaticRefineRegionLeftEdge"
            global StaticRefineRegionLeftEdge
            retval = []
            for i in range(MAX_STATIC_REGIONS):
                retval.append([])
                for j in range(MAX_DIMENSION):
                    retval[-1].append(StaticRefineRegionLeftEdge[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            global StaticRefineRegionLeftEdge
            for i in range(MAX_STATIC_REGIONS):
                for j in range(MAX_DIMENSION):
                    StaticRefineRegionLeftEdge[i][j] = val[i][j]

    property StaticRefineRegionRightEdge:
        def __get__(self):
            print "Returning a copy of StaticRefineRegionRightEdge"
            global StaticRefineRegionRightEdge
            retval = []
            for i in range(MAX_STATIC_REGIONS):
                retval.append([])
                for j in range(MAX_DIMENSION):
                    retval[-1].append(StaticRefineRegionRightEdge[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            global StaticRefineRegionRightEdge
            for i in range(MAX_STATIC_REGIONS):
                for j in range(MAX_DIMENSION):
                    StaticRefineRegionRightEdge[i][j] = val[i][j]

    property StaticRefineRegionLevel:
        def __get__(self):
            print "Returning a copy of StaticRefineRegionLevel"
            global StaticRefineRegionLevel
            retval = []
            for i in range(MAX_STATIC_REGIONS):
                retval.append(StaticRefineRegionLevel[i])
            return retval
        def __set__(self, val):
            cdef int i
            global StaticRefineRegionLevel
            for i in range(MAX_STATIC_REGIONS):
                StaticRefineRegionLevel[i] = val[i]

    property MyProcessorNumber:
        def __get__(self):
            global MyProcessorNumber
            return MyProcessorNumber
        def __set__(self, val):
            global MyProcessorNumber
            MyProcessorNumber = val

    property NumberOfProcessors:
        def __get__(self):
            global NumberOfProcessors
            return NumberOfProcessors
        def __set__(self, val):
            global NumberOfProcessors
            NumberOfProcessors = val

    property CommunicationTime:
        def __get__(self):
            global CommunicationTime
            return CommunicationTime
        def __set__(self, val):
            global CommunicationTime
            CommunicationTime = val

    property ParallelRootGridIO:
        def __get__(self):
            global ParallelRootGridIO
            return ParallelRootGridIO
        def __set__(self, val):
            global ParallelRootGridIO
            ParallelRootGridIO = val

    property ParallelParticleIO:
        def __get__(self):
            global ParallelParticleIO
            return ParallelParticleIO
        def __set__(self, val):
            global ParallelParticleIO
            ParallelParticleIO = val

    property Unigrid:
        def __get__(self):
            global Unigrid
            return Unigrid
        def __set__(self, val):
            global Unigrid
            Unigrid = val

    property CubeDumpEnabled:
        def __get__(self):
            global CubeDumpEnabled
            return CubeDumpEnabled
        def __set__(self, val):
            global CubeDumpEnabled
            CubeDumpEnabled = val

    property PartitionNestedGrids:
        def __get__(self):
            global PartitionNestedGrids
            return PartitionNestedGrids
        def __set__(self, val):
            global PartitionNestedGrids
            PartitionNestedGrids = val

    property ExtractFieldsOnly:
        def __get__(self):
            global ExtractFieldsOnly
            return ExtractFieldsOnly
        def __set__(self, val):
            global ExtractFieldsOnly
            ExtractFieldsOnly = val

    property First_Pass:
        def __get__(self):
            global First_Pass
            return First_Pass
        def __set__(self, val):
            global First_Pass
            First_Pass = val

    property UnigridTranspose:
        def __get__(self):
            global UnigridTranspose
            return UnigridTranspose
        def __set__(self, val):
            global UnigridTranspose
            UnigridTranspose = val

    property PythonSubcycleSkip:
        def __get__(self):
            global PythonSubcycleSkip
            return PythonSubcycleSkip
        def __set__(self, val):
            global PythonSubcycleSkip
            PythonSubcycleSkip = val

    property InlineHaloFinder:
        def __get__(self):
            global InlineHaloFinder
            return InlineHaloFinder
        def __set__(self, val):
            global InlineHaloFinder
            InlineHaloFinder = val

    property HaloFinderSubfind:
        def __get__(self):
            global HaloFinderSubfind
            return HaloFinderSubfind
        def __set__(self, val):
            global HaloFinderSubfind
            HaloFinderSubfind = val

    property HaloFinderOutputParticleList:
        def __get__(self):
            global HaloFinderOutputParticleList
            return HaloFinderOutputParticleList
        def __set__(self, val):
            global HaloFinderOutputParticleList
            HaloFinderOutputParticleList = val

    property HaloFinderMinimumSize:
        def __get__(self):
            global HaloFinderMinimumSize
            return HaloFinderMinimumSize
        def __set__(self, val):
            global HaloFinderMinimumSize
            HaloFinderMinimumSize = val

    property HaloFinderCycleSkip:
        def __get__(self):
            global HaloFinderCycleSkip
            return HaloFinderCycleSkip
        def __set__(self, val):
            global HaloFinderCycleSkip
            HaloFinderCycleSkip = val

    property HaloFinderLinkingLength:
        def __get__(self):
            global HaloFinderLinkingLength
            return HaloFinderLinkingLength
        def __set__(self, val):
            global HaloFinderLinkingLength
            HaloFinderLinkingLength = val

    property HaloFinderTimestep:
        def __get__(self):
            global HaloFinderTimestep
            return HaloFinderTimestep
        def __set__(self, val):
            global HaloFinderTimestep
            HaloFinderTimestep = val

    property HaloFinderLastTime:
        def __get__(self):
            global HaloFinderLastTime
            return HaloFinderLastTime
        def __set__(self, val):
            global HaloFinderLastTime
            HaloFinderLastTime = val

    property MinimumSlopeForRefinement:
        def __get__(self):
            print "Returning a copy of MinimumSlopeForRefinement"
            global MinimumSlopeForRefinement
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(MinimumSlopeForRefinement[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MinimumSlopeForRefinement
            for i in range(MAX_FLAGGING_METHODS):
                MinimumSlopeForRefinement[i] = val[i]

    property SlopeFlaggingFields:
        def __get__(self):
            print "Returning a copy of SlopeFlaggingFields"
            global SlopeFlaggingFields
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(SlopeFlaggingFields[i])
            return retval
        def __set__(self, val):
            cdef int i
            global SlopeFlaggingFields
            for i in range(MAX_FLAGGING_METHODS):
                SlopeFlaggingFields[i] = val[i]

    property MinimumOverDensityForRefinement:
        def __get__(self):
            print "Returning a copy of MinimumOverDensityForRefinement"
            global MinimumOverDensityForRefinement
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(MinimumOverDensityForRefinement[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MinimumOverDensityForRefinement
            for i in range(MAX_FLAGGING_METHODS):
                MinimumOverDensityForRefinement[i] = val[i]

    property MinimumMassForRefinement:
        def __get__(self):
            print "Returning a copy of MinimumMassForRefinement"
            global MinimumMassForRefinement
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(MinimumMassForRefinement[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MinimumMassForRefinement
            for i in range(MAX_FLAGGING_METHODS):
                MinimumMassForRefinement[i] = val[i]

    property MinimumMassForRefinementLevelExponent:
        def __get__(self):
            print "Returning a copy of MinimumMassForRefinementLevelExponent"
            global MinimumMassForRefinementLevelExponent
            retval = []
            for i in range(MAX_FLAGGING_METHODS):
                retval.append(MinimumMassForRefinementLevelExponent[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MinimumMassForRefinementLevelExponent
            for i in range(MAX_FLAGGING_METHODS):
                MinimumMassForRefinementLevelExponent[i] = val[i]

    property DepositPositionsParticleSmoothRadius:
        def __get__(self):
            global DepositPositionsParticleSmoothRadius
            return DepositPositionsParticleSmoothRadius
        def __set__(self, val):
            global DepositPositionsParticleSmoothRadius
            DepositPositionsParticleSmoothRadius = val

    property MinimumPressureJumpForRefinement:
        def __get__(self):
            global MinimumPressureJumpForRefinement
            return MinimumPressureJumpForRefinement
        def __set__(self, val):
            global MinimumPressureJumpForRefinement
            MinimumPressureJumpForRefinement = val

    property MinimumEnergyRatioForRefinement:
        def __get__(self):
            global MinimumEnergyRatioForRefinement
            return MinimumEnergyRatioForRefinement
        def __set__(self, val):
            global MinimumEnergyRatioForRefinement
            MinimumEnergyRatioForRefinement = val

    property RefineByJeansLengthSafetyFactor:
        def __get__(self):
            global RefineByJeansLengthSafetyFactor
            return RefineByJeansLengthSafetyFactor
        def __set__(self, val):
            global RefineByJeansLengthSafetyFactor
            RefineByJeansLengthSafetyFactor = val

    property MustRefineParticlesRefineToLevel:
        def __get__(self):
            global MustRefineParticlesRefineToLevel
            return MustRefineParticlesRefineToLevel
        def __set__(self, val):
            global MustRefineParticlesRefineToLevel
            MustRefineParticlesRefineToLevel = val

    property MinimumShearForRefinement:
        def __get__(self):
            global MinimumShearForRefinement
            return MinimumShearForRefinement
        def __set__(self, val):
            global MinimumShearForRefinement
            MinimumShearForRefinement = val

    property RefineByResistiveLengthSafetyFactor:
        def __get__(self):
            global RefineByResistiveLengthSafetyFactor
            return RefineByResistiveLengthSafetyFactor
        def __set__(self, val):
            global RefineByResistiveLengthSafetyFactor
            RefineByResistiveLengthSafetyFactor = val

    property NohProblemFullBox:
        def __get__(self):
            global NohProblemFullBox
            return NohProblemFullBox
        def __set__(self, val):
            global NohProblemFullBox
            NohProblemFullBox = val

    property ComovingCoordinates:
        def __get__(self):
            global ComovingCoordinates
            return ComovingCoordinates
        def __set__(self, val):
            global ComovingCoordinates
            ComovingCoordinates = val

    property StarParticleCreation:
        def __get__(self):
            global StarParticleCreation
            return StarParticleCreation
        def __set__(self, val):
            global StarParticleCreation
            StarParticleCreation = val

    property StarParticleFeedback:
        def __get__(self):
            global StarParticleFeedback
            return StarParticleFeedback
        def __set__(self, val):
            global StarParticleFeedback
            StarParticleFeedback = val

    property NumberOfParticleAttributes:
        def __get__(self):
            global NumberOfParticleAttributes
            return NumberOfParticleAttributes
        def __set__(self, val):
            global NumberOfParticleAttributes
            NumberOfParticleAttributes = val

    property AddParticleAttributes:
        def __get__(self):
            global AddParticleAttributes
            return AddParticleAttributes
        def __set__(self, val):
            global AddParticleAttributes
            AddParticleAttributes = val

    property TimeActionType:
        def __get__(self):
            print "Returning a copy of TimeActionType"
            global TimeActionType
            retval = []
            for i in range(MAX_TIME_ACTIONS):
                retval.append(TimeActionType[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TimeActionType
            for i in range(MAX_TIME_ACTIONS):
                TimeActionType[i] = val[i]

    property TimeActionTime:
        def __get__(self):
            print "Returning a copy of TimeActionTime"
            global TimeActionTime
            retval = []
            for i in range(MAX_TIME_ACTIONS):
                retval.append(TimeActionTime[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TimeActionTime
            for i in range(MAX_TIME_ACTIONS):
                TimeActionTime[i] = val[i]

    property TimeActionRedshift:
        def __get__(self):
            print "Returning a copy of TimeActionRedshift"
            global TimeActionRedshift
            retval = []
            for i in range(MAX_TIME_ACTIONS):
                retval.append(TimeActionRedshift[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TimeActionRedshift
            for i in range(MAX_TIME_ACTIONS):
                TimeActionRedshift[i] = val[i]

    property TimeActionParameter:
        def __get__(self):
            print "Returning a copy of TimeActionParameter"
            global TimeActionParameter
            retval = []
            for i in range(MAX_TIME_ACTIONS):
                retval.append(TimeActionParameter[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TimeActionParameter
            for i in range(MAX_TIME_ACTIONS):
                TimeActionParameter[i] = val[i]

    property TracerParticleOn:
        def __get__(self):
            global TracerParticleOn
            return TracerParticleOn
        def __set__(self, val):
            global TracerParticleOn
            TracerParticleOn = val

    property TracerParticleCreationSpacing:
        def __get__(self):
            global TracerParticleCreationSpacing
            return TracerParticleCreationSpacing
        def __set__(self, val):
            global TracerParticleCreationSpacing
            TracerParticleCreationSpacing = val

    property TracerParticleCreationLeftEdge:
        def __get__(self):
            print "Returning a copy of TracerParticleCreationLeftEdge"
            global TracerParticleCreationLeftEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(TracerParticleCreationLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TracerParticleCreationLeftEdge
            for i in range(MAX_DIMENSION):
                TracerParticleCreationLeftEdge[i] = val[i]

    property TracerParticleCreationRightEdge:
        def __get__(self):
            print "Returning a copy of TracerParticleCreationRightEdge"
            global TracerParticleCreationRightEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(TracerParticleCreationRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TracerParticleCreationRightEdge
            for i in range(MAX_DIMENSION):
                TracerParticleCreationRightEdge[i] = val[i]

    property ParticleTypeInFile:
        def __get__(self):
            global ParticleTypeInFile
            return ParticleTypeInFile
        def __set__(self, val):
            global ParticleTypeInFile
            ParticleTypeInFile = val

    property OutputParticleTypeGrouping:
        def __get__(self):
            global OutputParticleTypeGrouping
            return OutputParticleTypeGrouping
        def __set__(self, val):
            global OutputParticleTypeGrouping
            OutputParticleTypeGrouping = val

    property ExternalBoundaryIO:
        def __get__(self):
            global ExternalBoundaryIO
            return ExternalBoundaryIO
        def __set__(self, val):
            global ExternalBoundaryIO
            ExternalBoundaryIO = val

    property ExternalBoundaryTypeIO:
        def __get__(self):
            global ExternalBoundaryTypeIO
            return ExternalBoundaryTypeIO
        def __set__(self, val):
            global ExternalBoundaryTypeIO
            ExternalBoundaryTypeIO = val

    property ExternalBoundaryValueIO:
        def __get__(self):
            global ExternalBoundaryValueIO
            return ExternalBoundaryValueIO
        def __set__(self, val):
            global ExternalBoundaryValueIO
            ExternalBoundaryValueIO = val

    property ExternalBoundaryField:
        def __get__(self):
            global ExternalBoundaryField
            return ExternalBoundaryField
        def __set__(self, val):
            global ExternalBoundaryField
            ExternalBoundaryField = val

    property SimpleConstantBoundary:
        def __get__(self):
            global SimpleConstantBoundary
            return SimpleConstantBoundary
        def __set__(self, val):
            global SimpleConstantBoundary
            SimpleConstantBoundary = val

    property TaskMap:
        def __get__(self):
            print "Returning a copy of TaskMap"
            global TaskMap
            retval = []
            for i in range(MAX_NUMBER_OF_TASKS):
                retval.append(TaskMap[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TaskMap
            for i in range(MAX_NUMBER_OF_TASKS):
                TaskMap[i] = val[i]

    property NodeMem:
        def __get__(self):
            print "Returning a copy of NodeMem"
            global NodeMem
            retval = []
            for i in range(MAX_NUMBER_OF_NODES):
                retval.append(NodeMem[i])
            return retval
        def __set__(self, val):
            cdef int i
            global NodeMem
            for i in range(MAX_NUMBER_OF_NODES):
                NodeMem[i] = val[i]

    property NodeMap:
        def __get__(self):
            print "Returning a copy of NodeMap"
            global NodeMap
            retval = []
            for i in range(MAX_NUMBER_OF_NODES):
                retval.append(NodeMap[i])
            return retval
        def __set__(self, val):
            cdef int i
            global NodeMap
            for i in range(MAX_NUMBER_OF_NODES):
                NodeMap[i] = val[i]

    property LoadGridDataAtStart:
        def __get__(self):
            global LoadGridDataAtStart
            return LoadGridDataAtStart
        def __set__(self, val):
            global LoadGridDataAtStart
            LoadGridDataAtStart = val

    property GlobalCommunication:
        def __get__(self):
            global GlobalCommunication
            return GlobalCommunication
        def __set__(self, val):
            global GlobalCommunication
            GlobalCommunication = val

    property RecvComm:
        def __get__(self):
            global RecvComm
            return RecvComm
        def __set__(self, val):
            global RecvComm
            RecvComm = val

    property WaitComm:
        def __get__(self):
            global WaitComm
            return WaitComm
        def __set__(self, val):
            global WaitComm
            WaitComm = val

    property timer:
        def __get__(self):
            print "Returning a copy of timer"
            global timer
            retval = []
            for i in range(MAX_COUNTERS):
                retval.append(timer[i])
            return retval
        def __set__(self, val):
            cdef int i
            global timer
            for i in range(MAX_COUNTERS):
                timer[i] = val[i]

    property counter:
        def __get__(self):
            print "Returning a copy of counter"
            global counter
            retval = []
            for i in range(MAX_COUNTERS):
                retval.append(counter[i])
            return retval
        def __set__(self, val):
            cdef int i
            global counter
            for i in range(MAX_COUNTERS):
                counter[i] = val[i]

    property traceMEM:
        def __get__(self):
            global traceMEM
            return traceMEM
        def __set__(self, val):
            global traceMEM
            traceMEM = val

    property starttime:
        def __get__(self):
            global starttime
            return starttime
        def __set__(self, val):
            global starttime
            starttime = val

    property endtime:
        def __get__(self):
            global endtime
            return endtime
        def __set__(self, val):
            global endtime
            endtime = val

    property Start_Wall_Time:
        def __get__(self):
            global Start_Wall_Time
            return Start_Wall_Time
        def __set__(self, val):
            global Start_Wall_Time
            Start_Wall_Time = val

    property End_Wall_Time:
        def __get__(self):
            global End_Wall_Time
            return End_Wall_Time
        def __set__(self, val):
            global End_Wall_Time
            End_Wall_Time = val

    property WallTime:
        def __get__(self):
            global WallTime
            return WallTime
        def __set__(self, val):
            global WallTime
            WallTime = val

    property flagging_count:
        def __get__(self):
            global flagging_count
            return flagging_count
        def __set__(self, val):
            global flagging_count
            flagging_count = val

    property in_count:
        def __get__(self):
            global in_count
            return in_count
        def __set__(self, val):
            global in_count
            in_count = val

    property out_count:
        def __get__(self):
            global out_count
            return out_count
        def __set__(self, val):
            global out_count
            out_count = val

    property moving_count:
        def __get__(self):
            global moving_count
            return moving_count
        def __set__(self, val):
            global moving_count
            moving_count = val

    property flagging_pct:
        def __get__(self):
            global flagging_pct
            return flagging_pct
        def __set__(self, val):
            global flagging_pct
            flagging_pct = val

    property moving_pct:
        def __get__(self):
            global moving_pct
            return moving_pct
        def __set__(self, val):
            global moving_pct
            moving_pct = val

    property traceMPI:
        def __get__(self):
            global traceMPI
            return traceMPI
        def __set__(self, val):
            global traceMPI
            traceMPI = val

    property MovieDataField:
        def __get__(self):
            print "Returning a copy of MovieDataField"
            global MovieDataField
            retval = []
            for i in range(MAX_MOVIE_FIELDS):
                retval.append(MovieDataField[i])
            return retval
        def __set__(self, val):
            cdef int i
            global MovieDataField
            for i in range(MAX_MOVIE_FIELDS):
                MovieDataField[i] = val[i]

    property MovieSkipTimestep:
        def __get__(self):
            global MovieSkipTimestep
            return MovieSkipTimestep
        def __set__(self, val):
            global MovieSkipTimestep
            MovieSkipTimestep = val

    property Movie3DVolumes:
        def __get__(self):
            global Movie3DVolumes
            return Movie3DVolumes
        def __set__(self, val):
            global Movie3DVolumes
            Movie3DVolumes = val

    property MovieVertexCentered:
        def __get__(self):
            global MovieVertexCentered
            return MovieVertexCentered
        def __set__(self, val):
            global MovieVertexCentered
            MovieVertexCentered = val

    property NewMovieDumpNumber:
        def __get__(self):
            global NewMovieDumpNumber
            return NewMovieDumpNumber
        def __set__(self, val):
            global NewMovieDumpNumber
            NewMovieDumpNumber = val

    property NewMovieParticleOn:
        def __get__(self):
            global NewMovieParticleOn
            return NewMovieParticleOn
        def __set__(self, val):
            global NewMovieParticleOn
            NewMovieParticleOn = val

    property UseHydro:
        def __get__(self):
            global UseHydro
            return UseHydro
        def __set__(self, val):
            global UseHydro
            UseHydro = val

    property Coordinate:
        def __get__(self):
            global Coordinate
            return Coordinate
        def __set__(self, val):
            global Coordinate
            Coordinate = val

    property NSpecies:
        def __get__(self):
            global NSpecies
            return NSpecies
        def __set__(self, val):
            global NSpecies
            NSpecies = val

    property NColor:
        def __get__(self):
            global NColor
            return NColor
        def __set__(self, val):
            global NColor
            NColor = val

    property Theta_Limiter:
        def __get__(self):
            global Theta_Limiter
            return Theta_Limiter
        def __set__(self, val):
            global Theta_Limiter
            Theta_Limiter = val

    property RKOrder:
        def __get__(self):
            global RKOrder
            return RKOrder
        def __set__(self, val):
            global RKOrder
            RKOrder = val

    property UsePhysicalUnit:
        def __get__(self):
            global UsePhysicalUnit
            return UsePhysicalUnit
        def __set__(self, val):
            global UsePhysicalUnit
            UsePhysicalUnit = val

    property iden:
        def __get__(self):
            global iden
            return iden
        def __set__(self, val):
            global iden
            iden = val

    property ietot:
        def __get__(self):
            global ietot
            return ietot
        def __set__(self, val):
            global ietot
            ietot = val

    property ivx:
        def __get__(self):
            global ivx
            return ivx
        def __set__(self, val):
            global ivx
            ivx = val

    property ivy:
        def __get__(self):
            global ivy
            return ivy
        def __set__(self, val):
            global ivy
            ivy = val

    property ivz:
        def __get__(self):
            global ivz
            return ivz
        def __set__(self, val):
            global ivz
            ivz = val

    property iBx:
        def __get__(self):
            global iBx
            return iBx
        def __set__(self, val):
            global iBx
            iBx = val

    property iBy:
        def __get__(self):
            global iBy
            return iBy
        def __set__(self, val):
            global iBy
            iBy = val

    property iBz:
        def __get__(self):
            global iBz
            return iBz
        def __set__(self, val):
            global iBz
            iBz = val

    property iPhi:
        def __get__(self):
            global iPhi
            return iPhi
        def __set__(self, val):
            global iPhi
            iPhi = val

    property ieint:
        def __get__(self):
            global ieint
            return ieint
        def __set__(self, val):
            global ieint
            ieint = val

    property iD:
        def __get__(self):
            global iD
            return iD
        def __set__(self, val):
            global iD
            iD = val

    property iEtot:
        def __get__(self):
            global iEtot
            return iEtot
        def __set__(self, val):
            global iEtot
            iEtot = val

    property iS1:
        def __get__(self):
            global iS1
            return iS1
        def __set__(self, val):
            global iS1
            iS1 = val

    property iS2:
        def __get__(self):
            global iS2
            return iS2
        def __set__(self, val):
            global iS2
            iS2 = val

    property iS3:
        def __get__(self):
            global iS3
            return iS3
        def __set__(self, val):
            global iS3
            iS3 = val

    property iEint:
        def __get__(self):
            global iEint
            return iEint
        def __set__(self, val):
            global iEint
            iEint = val

    property SmallRho:
        def __get__(self):
            global SmallRho
            return SmallRho
        def __set__(self, val):
            global SmallRho
            SmallRho = val

    property SmallP:
        def __get__(self):
            global SmallP
            return SmallP
        def __set__(self, val):
            global SmallP
            SmallP = val

    property SmallEint:
        def __get__(self):
            global SmallEint
            return SmallEint
        def __set__(self, val):
            global SmallEint
            SmallEint = val

    property SmallT:
        def __get__(self):
            global SmallT
            return SmallT
        def __set__(self, val):
            global SmallT
            SmallT = val

    property MaximumAlvenSpeed:
        def __get__(self):
            global MaximumAlvenSpeed
            return MaximumAlvenSpeed
        def __set__(self, val):
            global MaximumAlvenSpeed
            MaximumAlvenSpeed = val

    property NEQ_HYDRO:
        def __get__(self):
            global NEQ_HYDRO
            return NEQ_HYDRO
        def __set__(self, val):
            global NEQ_HYDRO
            NEQ_HYDRO = val

    property NEQ_MHD:
        def __get__(self):
            global NEQ_MHD
            return NEQ_MHD
        def __set__(self, val):
            global NEQ_MHD
            NEQ_MHD = val

    property ReconstructionMethod:
        def __get__(self):
            global ReconstructionMethod
            return ReconstructionMethod
        def __set__(self, val):
            global ReconstructionMethod
            ReconstructionMethod = val

    property RiemannSolver:
        def __get__(self):
            global RiemannSolver
            return RiemannSolver
        def __set__(self, val):
            global RiemannSolver
            RiemannSolver = val

    property EOSType:
        def __get__(self):
            global EOSType
            return EOSType
        def __set__(self, val):
            global EOSType
            EOSType = val

    property EOSSoundSpeed:
        def __get__(self):
            global EOSSoundSpeed
            return EOSSoundSpeed
        def __set__(self, val):
            global EOSSoundSpeed
            EOSSoundSpeed = val

    property EOSCriticalDensity:
        def __get__(self):
            global EOSCriticalDensity
            return EOSCriticalDensity
        def __set__(self, val):
            global EOSCriticalDensity
            EOSCriticalDensity = val

    property EOSGamma:
        def __get__(self):
            global EOSGamma
            return EOSGamma
        def __set__(self, val):
            global EOSGamma
            EOSGamma = val

    property C_h:
        def __get__(self):
            global C_h
            return C_h
        def __set__(self, val):
            global C_h
            C_h = val

    property C_p:
        def __get__(self):
            global C_p
            return C_p
        def __set__(self, val):
            global C_p
            C_p = val

    property UseConstantAcceleration:
        def __get__(self):
            global UseConstantAcceleration
            return UseConstantAcceleration
        def __set__(self, val):
            global UseConstantAcceleration
            UseConstantAcceleration = val

    property ConstantAcceleration:
        def __get__(self):
            print "Returning a copy of ConstantAcceleration"
            global ConstantAcceleration
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(ConstantAcceleration[i])
            return retval
        def __set__(self, val):
            cdef int i
            global ConstantAcceleration
            for i in range(MAX_DIMENSION):
                ConstantAcceleration[i] = val[i]

    property Mu:
        def __get__(self):
            global Mu
            return Mu
        def __set__(self, val):
            global Mu
            Mu = val

    property ExternalGravity:
        def __get__(self):
            global ExternalGravity
            return ExternalGravity
        def __set__(self, val):
            global ExternalGravity
            ExternalGravity = val

    property StringKick:
        def __get__(self):
            global StringKick
            return StringKick
        def __set__(self, val):
            global StringKick
            StringKick = val

    property UseFloor:
        def __get__(self):
            global UseFloor
            return UseFloor
        def __set__(self, val):
            global UseFloor
            UseFloor = val

    property UseViscosity:
        def __get__(self):
            global UseViscosity
            return UseViscosity
        def __set__(self, val):
            global UseViscosity
            UseViscosity = val

    property UseAmbipolarDiffusion:
        def __get__(self):
            global UseAmbipolarDiffusion
            return UseAmbipolarDiffusion
        def __set__(self, val):
            global UseAmbipolarDiffusion
            UseAmbipolarDiffusion = val

    property UseResistivity:
        def __get__(self):
            global UseResistivity
            return UseResistivity
        def __set__(self, val):
            global UseResistivity
            UseResistivity = val

    property UseH2OnDust:
        def __get__(self):
            global UseH2OnDust
            return UseH2OnDust
        def __set__(self, val):
            global UseH2OnDust
            UseH2OnDust = val

    property PhotoelectricHeating:
        def __get__(self):
            global PhotoelectricHeating
            return PhotoelectricHeating
        def __set__(self, val):
            global PhotoelectricHeating
            PhotoelectricHeating = val

    property CoolingCutOffDensity1:
        def __get__(self):
            global CoolingCutOffDensity1
            return CoolingCutOffDensity1
        def __set__(self, val):
            global CoolingCutOffDensity1
            CoolingCutOffDensity1 = val

    property CoolingCutOffDensity2:
        def __get__(self):
            global CoolingCutOffDensity2
            return CoolingCutOffDensity2
        def __set__(self, val):
            global CoolingCutOffDensity2
            CoolingCutOffDensity2 = val

    property CoolingPowerCutOffDensity1:
        def __get__(self):
            global CoolingPowerCutOffDensity1
            return CoolingPowerCutOffDensity1
        def __set__(self, val):
            global CoolingPowerCutOffDensity1
            CoolingPowerCutOffDensity1 = val

    property CoolingPowerCutOffDensity2:
        def __get__(self):
            global CoolingPowerCutOffDensity2
            return CoolingPowerCutOffDensity2
        def __set__(self, val):
            global CoolingPowerCutOffDensity2
            CoolingPowerCutOffDensity2 = val

    property CoolingCutOffTemperature:
        def __get__(self):
            global CoolingCutOffTemperature
            return CoolingCutOffTemperature
        def __set__(self, val):
            global CoolingCutOffTemperature
            CoolingCutOffTemperature = val

    property CoolingModel:
        def __get__(self):
            global CoolingModel
            return CoolingModel
        def __set__(self, val):
            global CoolingModel
            CoolingModel = val

    property HaloMass:
        def __get__(self):
            global HaloMass
            return HaloMass
        def __set__(self, val):
            global HaloMass
            HaloMass = val

    property HaloConcentration:
        def __get__(self):
            global HaloConcentration
            return HaloConcentration
        def __set__(self, val):
            global HaloConcentration
            HaloConcentration = val

    property HaloRedshift:
        def __get__(self):
            global HaloRedshift
            return HaloRedshift
        def __set__(self, val):
            global HaloRedshift
            HaloRedshift = val

    property HaloCentralDensity:
        def __get__(self):
            global HaloCentralDensity
            return HaloCentralDensity
        def __set__(self, val):
            global HaloCentralDensity
            HaloCentralDensity = val

    property HaloVirialRadius:
        def __get__(self):
            global HaloVirialRadius
            return HaloVirialRadius
        def __set__(self, val):
            global HaloVirialRadius
            HaloVirialRadius = val

    property ExternalGravityDensity:
        def __get__(self):
            global ExternalGravityDensity
            return ExternalGravityDensity
        def __set__(self, val):
            global ExternalGravityDensity
            ExternalGravityDensity = val

    property ExternalGravityRadius:
        def __get__(self):
            global ExternalGravityRadius
            return ExternalGravityRadius
        def __set__(self, val):
            global ExternalGravityRadius
            ExternalGravityRadius = val

    property UseDivergenceCleaning:
        def __get__(self):
            global UseDivergenceCleaning
            return UseDivergenceCleaning
        def __set__(self, val):
            global UseDivergenceCleaning
            UseDivergenceCleaning = val

    property DivergenceCleaningBoundaryBuffer:
        def __get__(self):
            global DivergenceCleaningBoundaryBuffer
            return DivergenceCleaningBoundaryBuffer
        def __set__(self, val):
            global DivergenceCleaningBoundaryBuffer
            DivergenceCleaningBoundaryBuffer = val

    property DivergenceCleaningThreshold:
        def __get__(self):
            global DivergenceCleaningThreshold
            return DivergenceCleaningThreshold
        def __set__(self, val):
            global DivergenceCleaningThreshold
            DivergenceCleaningThreshold = val

    property PoissonApproximationThreshold:
        def __get__(self):
            global PoissonApproximationThreshold
            return PoissonApproximationThreshold
        def __set__(self, val):
            global PoissonApproximationThreshold
            PoissonApproximationThreshold = val

    property ShiningParticleID:
        def __get__(self):
            global ShiningParticleID
            return ShiningParticleID
        def __set__(self, val):
            global ShiningParticleID
            ShiningParticleID = val

    property SinkMergeDistance:
        def __get__(self):
            global SinkMergeDistance
            return SinkMergeDistance
        def __set__(self, val):
            global SinkMergeDistance
            SinkMergeDistance = val

    property SinkMergeMass:
        def __get__(self):
            global SinkMergeMass
            return SinkMergeMass
        def __set__(self, val):
            global SinkMergeMass
            SinkMergeMass = val

    property TotalSinkMass:
        def __get__(self):
            global TotalSinkMass
            return TotalSinkMass
        def __set__(self, val):
            global TotalSinkMass
            TotalSinkMass = val

    property StellarWindFeedback:
        def __get__(self):
            global StellarWindFeedback
            return StellarWindFeedback
        def __set__(self, val):
            global StellarWindFeedback
            StellarWindFeedback = val

    property StellarWindTurnOnMass:
        def __get__(self):
            global StellarWindTurnOnMass
            return StellarWindTurnOnMass
        def __set__(self, val):
            global StellarWindTurnOnMass
            StellarWindTurnOnMass = val

    property NBodyDirectSummation:
        def __get__(self):
            global NBodyDirectSummation
            return NBodyDirectSummation
        def __set__(self, val):
            global NBodyDirectSummation
            NBodyDirectSummation = val

    property UseDrivingField:
        def __get__(self):
            global UseDrivingField
            return UseDrivingField
        def __set__(self, val):
            global UseDrivingField
            UseDrivingField = val

    property DrivingEfficiency:
        def __get__(self):
            global DrivingEfficiency
            return DrivingEfficiency
        def __set__(self, val):
            global DrivingEfficiency
            DrivingEfficiency = val

    property UseCUDA:
        def __get__(self):
            global UseCUDA
            return UseCUDA
        def __set__(self, val):
            global UseCUDA
            UseCUDA = val

    property ran1_init:
        def __get__(self):
            global ran1_init
            return ran1_init
        def __set__(self, val):
            global ran1_init
            ran1_init = val

    property MemoryLimit:
        def __get__(self):
            global MemoryLimit
            return MemoryLimit
        def __set__(self, val):
            global MemoryLimit
            MemoryLimit = val

    property MetalCooling:
        def __get__(self):
            global MetalCooling
            return MetalCooling
        def __set__(self, val):
            global MetalCooling
            MetalCooling = val

    property CIECooling:
        def __get__(self):
            global CIECooling
            return CIECooling
        def __set__(self, val):
            global CIECooling
            CIECooling = val

    property H2OpticalDepthApproximation:
        def __get__(self):
            global H2OpticalDepthApproximation
            return H2OpticalDepthApproximation
        def __set__(self, val):
            global H2OpticalDepthApproximation
            H2OpticalDepthApproximation = val

    property RadiativeTransfer:
        def __get__(self):
            global RadiativeTransfer
            return RadiativeTransfer
        def __set__(self, val):
            global RadiativeTransfer
            RadiativeTransfer = val

    property PhotonTime:
        def __get__(self):
            global PhotonTime
            return PhotonTime
        def __set__(self, val):
            global PhotonTime
            PhotonTime = val

    property dtPhoton:
        def __get__(self):
            global dtPhoton
            return dtPhoton
        def __set__(self, val):
            global dtPhoton
            dtPhoton = val

    property RadiativeTransferSourceRadius:
        def __get__(self):
            global RadiativeTransferSourceRadius
            return RadiativeTransferSourceRadius
        def __set__(self, val):
            global RadiativeTransferSourceRadius
            RadiativeTransferSourceRadius = val

    property RadiativeTransferPropagationSpeedFraction:
        def __get__(self):
            global RadiativeTransferPropagationSpeedFraction
            return RadiativeTransferPropagationSpeedFraction
        def __set__(self, val):
            global RadiativeTransferPropagationSpeedFraction
            RadiativeTransferPropagationSpeedFraction = val

    property RadiativeTransferPropagationDistance:
        def __get__(self):
            global RadiativeTransferPropagationDistance
            return RadiativeTransferPropagationDistance
        def __set__(self, val):
            global RadiativeTransferPropagationDistance
            RadiativeTransferPropagationDistance = val

    property RadiativeTransferSplitPhotonRadius:
        def __get__(self):
            global RadiativeTransferSplitPhotonRadius
            return RadiativeTransferSplitPhotonRadius
        def __set__(self, val):
            global RadiativeTransferSplitPhotonRadius
            RadiativeTransferSplitPhotonRadius = val

    property RadiativeTransferRaysPerCell:
        def __get__(self):
            global RadiativeTransferRaysPerCell
            return RadiativeTransferRaysPerCell
        def __set__(self, val):
            global RadiativeTransferRaysPerCell
            RadiativeTransferRaysPerCell = val

    property RadiativeTransferInitialHEALPixLevel:
        def __get__(self):
            global RadiativeTransferInitialHEALPixLevel
            return RadiativeTransferInitialHEALPixLevel
        def __set__(self, val):
            global RadiativeTransferInitialHEALPixLevel
            RadiativeTransferInitialHEALPixLevel = val

    property RadiativeTransferPhotonEscapeRadius:
        def __get__(self):
            global RadiativeTransferPhotonEscapeRadius
            return RadiativeTransferPhotonEscapeRadius
        def __set__(self, val):
            global RadiativeTransferPhotonEscapeRadius
            RadiativeTransferPhotonEscapeRadius = val

    property RadiativeTransferInterpolateField:
        def __get__(self):
            global RadiativeTransferInterpolateField
            return RadiativeTransferInterpolateField
        def __set__(self, val):
            global RadiativeTransferInterpolateField
            RadiativeTransferInterpolateField = val

    property RadiativeTransferTimestepVelocityLimit:
        def __get__(self):
            global RadiativeTransferTimestepVelocityLimit
            return RadiativeTransferTimestepVelocityLimit
        def __set__(self, val):
            global RadiativeTransferTimestepVelocityLimit
            RadiativeTransferTimestepVelocityLimit = val

    property RadiativeTransferSourceClustering:
        def __get__(self):
            global RadiativeTransferSourceClustering
            return RadiativeTransferSourceClustering
        def __set__(self, val):
            global RadiativeTransferSourceClustering
            RadiativeTransferSourceClustering = val

    property RadiativeTransferPhotonMergeRadius:
        def __get__(self):
            global RadiativeTransferPhotonMergeRadius
            return RadiativeTransferPhotonMergeRadius
        def __set__(self, val):
            global RadiativeTransferPhotonMergeRadius
            RadiativeTransferPhotonMergeRadius = val

    property RadiationPressure:
        def __get__(self):
            global RadiationPressure
            return RadiationPressure
        def __set__(self, val):
            global RadiationPressure
            RadiationPressure = val

    property RadiativeTransferOpticallyThinH2:
        def __get__(self):
            global RadiativeTransferOpticallyThinH2
            return RadiativeTransferOpticallyThinH2
        def __set__(self, val):
            global RadiativeTransferOpticallyThinH2
            RadiativeTransferOpticallyThinH2 = val

    property RadiativeTransferPeriodicBoundary:
        def __get__(self):
            global RadiativeTransferPeriodicBoundary
            return RadiativeTransferPeriodicBoundary
        def __set__(self, val):
            global RadiativeTransferPeriodicBoundary
            RadiativeTransferPeriodicBoundary = val

    property EscapedPhotonCount:
        def __get__(self):
            print "Returning a copy of EscapedPhotonCount"
            global EscapedPhotonCount
            retval = []
            for i in range(4):
                retval.append(EscapedPhotonCount[i])
            return retval
        def __set__(self, val):
            cdef int i
            global EscapedPhotonCount
            for i in range(4):
                EscapedPhotonCount[i] = val[i]

    property TotalEscapedPhotonCount:
        def __get__(self):
            print "Returning a copy of TotalEscapedPhotonCount"
            global TotalEscapedPhotonCount
            retval = []
            for i in range(4):
                retval.append(TotalEscapedPhotonCount[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TotalEscapedPhotonCount
            for i in range(4):
                TotalEscapedPhotonCount[i] = val[i]

    property FieldsToInterpolate:
        def __get__(self):
            print "Returning a copy of FieldsToInterpolate"
            global FieldsToInterpolate
            retval = []
            for i in range(MAX_NUMBER_OF_BARYON_FIELDS):
                retval.append(FieldsToInterpolate[i])
            return retval
        def __set__(self, val):
            cdef int i
            global FieldsToInterpolate
            for i in range(MAX_NUMBER_OF_BARYON_FIELDS):
                FieldsToInterpolate[i] = val[i]

    property RadiativeTransferCoupledRateSolver:
        def __get__(self):
            global RadiativeTransferCoupledRateSolver
            return RadiativeTransferCoupledRateSolver
        def __set__(self, val):
            global RadiativeTransferCoupledRateSolver
            RadiativeTransferCoupledRateSolver = val

    property AngularVelocity:
        def __get__(self):
            global AngularVelocity
            return AngularVelocity
        def __set__(self, val):
            global AngularVelocity
            AngularVelocity = val

    property VelocityGradient:
        def __get__(self):
            global VelocityGradient
            return VelocityGradient
        def __set__(self, val):
            global VelocityGradient
            VelocityGradient = val

    property ShearingBoundaryDirection:
        def __get__(self):
            global ShearingBoundaryDirection
            return ShearingBoundaryDirection
        def __set__(self, val):
            global ShearingBoundaryDirection
            ShearingBoundaryDirection = val

    property ShearingVelocityDirection:
        def __get__(self):
            global ShearingVelocityDirection
            return ShearingVelocityDirection
        def __set__(self, val):
            global ShearingVelocityDirection
            ShearingVelocityDirection = val

    property ShearingOtherDirection:
        def __get__(self):
            global ShearingOtherDirection
            return ShearingOtherDirection
        def __set__(self, val):
            global ShearingOtherDirection
            ShearingOtherDirection = val

    property useMHD:
        def __get__(self):
            global useMHD
            return useMHD
        def __set__(self, val):
            global useMHD
            useMHD = val

    property TopGridDx:
        def __get__(self):
            print "Returning a copy of TopGridDx"
            global TopGridDx
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(TopGridDx[i])
            return retval
        def __set__(self, val):
            cdef int i
            global TopGridDx
            for i in range(MAX_DIMENSION):
                TopGridDx[i] = val[i]

    property ShearingBoxProblemType:
        def __get__(self):
            global ShearingBoxProblemType
            return ShearingBoxProblemType
        def __set__(self, val):
            global ShearingBoxProblemType
            ShearingBoxProblemType = val

    property IsothermalSoundSpeed:
        def __get__(self):
            global IsothermalSoundSpeed
            return IsothermalSoundSpeed
        def __set__(self, val):
            global IsothermalSoundSpeed
            IsothermalSoundSpeed = val

    property RefineByJeansLengthUnits:
        def __get__(self):
            global RefineByJeansLengthUnits
            return RefineByJeansLengthUnits
        def __set__(self, val):
            global RefineByJeansLengthUnits
            RefineByJeansLengthUnits = val

    property GlobalLengthUnits:
        def __get__(self):
            global GlobalLengthUnits
            return GlobalLengthUnits
        def __set__(self, val):
            global GlobalLengthUnits
            GlobalLengthUnits = val

    property GlobalMassUnits:
        def __get__(self):
            global GlobalMassUnits
            return GlobalMassUnits
        def __set__(self, val):
            global GlobalMassUnits
            GlobalMassUnits = val

    property GlobalDensityUnits:
        def __get__(self):
            global GlobalDensityUnits
            return GlobalDensityUnits
        def __set__(self, val):
            global GlobalDensityUnits
            GlobalDensityUnits = val

    property GlobalTimeUnits:
        def __get__(self):
            global GlobalTimeUnits
            return GlobalTimeUnits
        def __set__(self, val):
            global GlobalTimeUnits
            GlobalTimeUnits = val

    property flow_trace_on:
        def __get__(self):
            global flow_trace_on
            return flow_trace_on
        def __set__(self, val):
            global flow_trace_on
            flow_trace_on = val

    property flow_trace_level:
        def __get__(self):
            global flow_trace_level
            return flow_trace_level
        def __set__(self, val):
            global flow_trace_level
            flow_trace_level = val

    property CommunicationDirection:
        def __get__(self):
            global CommunicationDirection
            return CommunicationDirection
        def __set__(self, val):
            global CommunicationDirection
            CommunicationDirection = val

    property HubbleConstantNow:
        def __get__(self):
            global HubbleConstantNow
            return HubbleConstantNow
        def __set__(self, val):
            global HubbleConstantNow
            HubbleConstantNow = val

    property OmegaMatterNow:
        def __get__(self):
            global OmegaMatterNow
            return OmegaMatterNow
        def __set__(self, val):
            global OmegaMatterNow
            OmegaMatterNow = val

    property OmegaLambdaNow:
        def __get__(self):
            global OmegaLambdaNow
            return OmegaLambdaNow
        def __set__(self, val):
            global OmegaLambdaNow
            OmegaLambdaNow = val

    property ComovingBoxSize:
        def __get__(self):
            global ComovingBoxSize
            return ComovingBoxSize
        def __set__(self, val):
            global ComovingBoxSize
            ComovingBoxSize = val

    property MaxExpansionRate:
        def __get__(self):
            global MaxExpansionRate
            return MaxExpansionRate
        def __set__(self, val):
            global MaxExpansionRate
            MaxExpansionRate = val

    property InitialTimeInCodeUnits:
        def __get__(self):
            global InitialTimeInCodeUnits
            return InitialTimeInCodeUnits
        def __set__(self, val):
            global InitialTimeInCodeUnits
            InitialTimeInCodeUnits = val

    property InitialRedshift:
        def __get__(self):
            global InitialRedshift
            return InitialRedshift
        def __set__(self, val):
            global InitialRedshift
            InitialRedshift = val

    property CosmologyOutputRedshift:
        def __get__(self):
            print "Returning a copy of CosmologyOutputRedshift"
            global CosmologyOutputRedshift
            retval = []
            for i in range(MAX_NUMBER_OF_OUTPUT_REDSHIFTS):
                retval.append(CosmologyOutputRedshift[i])
            return retval
        def __set__(self, val):
            cdef int i
            global CosmologyOutputRedshift
            for i in range(MAX_NUMBER_OF_OUTPUT_REDSHIFTS):
                CosmologyOutputRedshift[i] = val[i]

    property CosmologyOutputRedshiftTime:
        def __get__(self):
            print "Returning a copy of CosmologyOutputRedshiftTime"
            global CosmologyOutputRedshiftTime
            retval = []
            for i in range(MAX_NUMBER_OF_OUTPUT_REDSHIFTS):
                retval.append(CosmologyOutputRedshiftTime[i])
            return retval
        def __set__(self, val):
            cdef int i
            global CosmologyOutputRedshiftTime
            for i in range(MAX_NUMBER_OF_OUTPUT_REDSHIFTS):
                CosmologyOutputRedshiftTime[i] = val[i]

    property NumberOfStarParticles:
        def __get__(self):
            global NumberOfStarParticles
            return NumberOfStarParticles
        def __set__(self, val):
            global NumberOfStarParticles
            NumberOfStarParticles = val

    property TotalNumberOfStars:
        def __get__(self):
            global TotalNumberOfStars
            return TotalNumberOfStars
        def __set__(self, val):
            global TotalNumberOfStars
            TotalNumberOfStars = val

    property StarMakerOverDensityThreshold:
        def __get__(self):
            global StarMakerOverDensityThreshold
            return StarMakerOverDensityThreshold
        def __set__(self, val):
            global StarMakerOverDensityThreshold
            StarMakerOverDensityThreshold = val

    property StarMakerMassEfficiency:
        def __get__(self):
            global StarMakerMassEfficiency
            return StarMakerMassEfficiency
        def __set__(self, val):
            global StarMakerMassEfficiency
            StarMakerMassEfficiency = val

    property StarMakerMinimumMass:
        def __get__(self):
            global StarMakerMinimumMass
            return StarMakerMinimumMass
        def __set__(self, val):
            global StarMakerMinimumMass
            StarMakerMinimumMass = val

    property StarMakerMinimumDynamicalTime:
        def __get__(self):
            global StarMakerMinimumDynamicalTime
            return StarMakerMinimumDynamicalTime
        def __set__(self, val):
            global StarMakerMinimumDynamicalTime
            StarMakerMinimumDynamicalTime = val

    property StarMassEjectionFraction:
        def __get__(self):
            global StarMassEjectionFraction
            return StarMassEjectionFraction
        def __set__(self, val):
            global StarMassEjectionFraction
            StarMassEjectionFraction = val

    property StarMetalYield:
        def __get__(self):
            global StarMetalYield
            return StarMetalYield
        def __set__(self, val):
            global StarMetalYield
            StarMetalYield = val

    property StarEnergyToThermalFeedback:
        def __get__(self):
            global StarEnergyToThermalFeedback
            return StarEnergyToThermalFeedback
        def __set__(self, val):
            global StarEnergyToThermalFeedback
            StarEnergyToThermalFeedback = val

    property StarEnergyToStellarUV:
        def __get__(self):
            global StarEnergyToStellarUV
            return StarEnergyToStellarUV
        def __set__(self, val):
            global StarEnergyToStellarUV
            StarEnergyToStellarUV = val

    property StarEnergyToQuasarUV:
        def __get__(self):
            global StarEnergyToQuasarUV
            return StarEnergyToQuasarUV
        def __set__(self, val):
            global StarEnergyToQuasarUV
            StarEnergyToQuasarUV = val

    property PopIIIStarMass:
        def __get__(self):
            global PopIIIStarMass
            return PopIIIStarMass
        def __set__(self, val):
            global PopIIIStarMass
            PopIIIStarMass = val

    property PopIIIBlackHoles:
        def __get__(self):
            global PopIIIBlackHoles
            return PopIIIBlackHoles
        def __set__(self, val):
            global PopIIIBlackHoles
            PopIIIBlackHoles = val

    property PopIIIBHLuminosityEfficiency:
        def __get__(self):
            global PopIIIBHLuminosityEfficiency
            return PopIIIBHLuminosityEfficiency
        def __set__(self, val):
            global PopIIIBHLuminosityEfficiency
            PopIIIBHLuminosityEfficiency = val

    property PopIIIOverDensityThreshold:
        def __get__(self):
            global PopIIIOverDensityThreshold
            return PopIIIOverDensityThreshold
        def __set__(self, val):
            global PopIIIOverDensityThreshold
            PopIIIOverDensityThreshold = val

    property PopIIIH2CriticalFraction:
        def __get__(self):
            global PopIIIH2CriticalFraction
            return PopIIIH2CriticalFraction
        def __set__(self, val):
            global PopIIIH2CriticalFraction
            PopIIIH2CriticalFraction = val

    property PopIIIMetalCriticalFraction:
        def __get__(self):
            global PopIIIMetalCriticalFraction
            return PopIIIMetalCriticalFraction
        def __set__(self, val):
            global PopIIIMetalCriticalFraction
            PopIIIMetalCriticalFraction = val

    property PopIIISupernovaRadius:
        def __get__(self):
            global PopIIISupernovaRadius
            return PopIIISupernovaRadius
        def __set__(self, val):
            global PopIIISupernovaRadius
            PopIIISupernovaRadius = val

    property PopIIISupernovaUseColour:
        def __get__(self):
            global PopIIISupernovaUseColour
            return PopIIISupernovaUseColour
        def __set__(self, val):
            global PopIIISupernovaUseColour
            PopIIISupernovaUseColour = val

    property StarClusterUseMetalField:
        def __get__(self):
            global StarClusterUseMetalField
            return StarClusterUseMetalField
        def __set__(self, val):
            global StarClusterUseMetalField
            StarClusterUseMetalField = val

    property StarClusterMinDynamicalTime:
        def __get__(self):
            global StarClusterMinDynamicalTime
            return StarClusterMinDynamicalTime
        def __set__(self, val):
            global StarClusterMinDynamicalTime
            StarClusterMinDynamicalTime = val

    property StarClusterIonizingLuminosity:
        def __get__(self):
            global StarClusterIonizingLuminosity
            return StarClusterIonizingLuminosity
        def __set__(self, val):
            global StarClusterIonizingLuminosity
            StarClusterIonizingLuminosity = val

    property StarClusterSNEnergy:
        def __get__(self):
            global StarClusterSNEnergy
            return StarClusterSNEnergy
        def __set__(self, val):
            global StarClusterSNEnergy
            StarClusterSNEnergy = val

    property StarClusterSNRadius:
        def __get__(self):
            global StarClusterSNRadius
            return StarClusterSNRadius
        def __set__(self, val):
            global StarClusterSNRadius
            StarClusterSNRadius = val

    property StarClusterFormEfficiency:
        def __get__(self):
            global StarClusterFormEfficiency
            return StarClusterFormEfficiency
        def __set__(self, val):
            global StarClusterFormEfficiency
            StarClusterFormEfficiency = val

    property StarClusterMinimumMass:
        def __get__(self):
            global StarClusterMinimumMass
            return StarClusterMinimumMass
        def __set__(self, val):
            global StarClusterMinimumMass
            StarClusterMinimumMass = val

    property StarClusterCombineRadius:
        def __get__(self):
            global StarClusterCombineRadius
            return StarClusterCombineRadius
        def __set__(self, val):
            global StarClusterCombineRadius
            StarClusterCombineRadius = val

    property StarClusterRegionLeftEdge:
        def __get__(self):
            print "Returning a copy of StarClusterRegionLeftEdge"
            global StarClusterRegionLeftEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(StarClusterRegionLeftEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global StarClusterRegionLeftEdge
            for i in range(MAX_DIMENSION):
                StarClusterRegionLeftEdge[i] = val[i]

    property StarClusterRegionRightEdge:
        def __get__(self):
            print "Returning a copy of StarClusterRegionRightEdge"
            global StarClusterRegionRightEdge
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append(StarClusterRegionRightEdge[i])
            return retval
        def __set__(self, val):
            cdef int i
            global StarClusterRegionRightEdge
            for i in range(MAX_DIMENSION):
                StarClusterRegionRightEdge[i] = val[i]

    property MBHMinDynamicalTime:
        def __get__(self):
            global MBHMinDynamicalTime
            return MBHMinDynamicalTime
        def __set__(self, val):
            global MBHMinDynamicalTime
            MBHMinDynamicalTime = val

    property MBHMinimumMass:
        def __get__(self):
            global MBHMinimumMass
            return MBHMinimumMass
        def __set__(self, val):
            global MBHMinimumMass
            MBHMinimumMass = val

    property MBHFeedbackThermal:
        def __get__(self):
            global MBHFeedbackThermal
            return MBHFeedbackThermal
        def __set__(self, val):
            global MBHFeedbackThermal
            MBHFeedbackThermal = val

    property MBHFeedbackRadius:
        def __get__(self):
            global MBHFeedbackRadius
            return MBHFeedbackRadius
        def __set__(self, val):
            global MBHFeedbackRadius
            MBHFeedbackRadius = val

    property MBHFeedbackRadiativeEfficiency:
        def __get__(self):
            global MBHFeedbackRadiativeEfficiency
            return MBHFeedbackRadiativeEfficiency
        def __set__(self, val):
            global MBHFeedbackRadiativeEfficiency
            MBHFeedbackRadiativeEfficiency = val

    property MBHFeedbackThermalCoupling:
        def __get__(self):
            global MBHFeedbackThermalCoupling
            return MBHFeedbackThermalCoupling
        def __set__(self, val):
            global MBHFeedbackThermalCoupling
            MBHFeedbackThermalCoupling = val

    property MBHCombineRadius:
        def __get__(self):
            global MBHCombineRadius
            return MBHCombineRadius
        def __set__(self, val):
            global MBHCombineRadius
            MBHCombineRadius = val

    property MBHIonizingLuminosity:
        def __get__(self):
            global MBHIonizingLuminosity
            return MBHIonizingLuminosity
        def __set__(self, val):
            global MBHIonizingLuminosity
            MBHIonizingLuminosity = val

    property minStarLifetime:
        def __get__(self):
            global minStarLifetime
            return minStarLifetime
        def __set__(self, val):
            global minStarLifetime
            minStarLifetime = val

    property LastSupernovaTime:
        def __get__(self):
            global LastSupernovaTime
            return LastSupernovaTime
        def __set__(self, val):
            global LastSupernovaTime
            LastSupernovaTime = val

    property CommunicationReceiveCurrentDependsOn:
        def __get__(self):
            global CommunicationReceiveCurrentDependsOn
            return CommunicationReceiveCurrentDependsOn
        def __set__(self, val):
            global CommunicationReceiveCurrentDependsOn
            CommunicationReceiveCurrentDependsOn = val

    property CommunicationReceiveIndex:
        def __get__(self):
            global CommunicationReceiveIndex
            return CommunicationReceiveIndex
        def __set__(self, val):
            global CommunicationReceiveIndex
            CommunicationReceiveIndex = val

    property CommunicationReceiveCallType:
        def __get__(self):
            print "Returning a copy of CommunicationReceiveCallType"
            global CommunicationReceiveCallType
            retval = []
            for i in range(MAX_RECEIVE_BUFFERS):
                retval.append(CommunicationReceiveCallType[i])
            return retval
        def __set__(self, val):
            cdef int i
            global CommunicationReceiveCallType
            for i in range(MAX_RECEIVE_BUFFERS):
                CommunicationReceiveCallType[i] = val[i]

    property CommunicationReceiveDependsOn:
        def __get__(self):
            print "Returning a copy of CommunicationReceiveDependsOn"
            global CommunicationReceiveDependsOn
            retval = []
            for i in range(MAX_RECEIVE_BUFFERS):
                retval.append(CommunicationReceiveDependsOn[i])
            return retval
        def __set__(self, val):
            cdef int i
            global CommunicationReceiveDependsOn
            for i in range(MAX_RECEIVE_BUFFERS):
                CommunicationReceiveDependsOn[i] = val[i]

    property CommunicationReceiveArgument:
        def __get__(self):
            print "Returning a copy of CommunicationReceiveArgument"
            global CommunicationReceiveArgument
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append([])
                for j in range(MAX_RECEIVE_BUFFERS):
                    retval[-1].append(CommunicationReceiveArgument[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            global CommunicationReceiveArgument
            for i in range(MAX_DIMENSION):
                for j in range(MAX_RECEIVE_BUFFERS):
                    CommunicationReceiveArgument[i][j] = val[i][j]

    property CommunicationReceiveArgumentInt:
        def __get__(self):
            print "Returning a copy of CommunicationReceiveArgumentInt"
            global CommunicationReceiveArgumentInt
            retval = []
            for i in range(MAX_DIMENSION):
                retval.append([])
                for j in range(MAX_RECEIVE_BUFFERS):
                    retval[-1].append(CommunicationReceiveArgumentInt[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            global CommunicationReceiveArgumentInt
            for i in range(MAX_DIMENSION):
                for j in range(MAX_RECEIVE_BUFFERS):
                    CommunicationReceiveArgumentInt[i][j] = val[i][j]

# END GLOBAL DATA HANDLER
