from mpi4py import MPI
import enzo_wrap as ew
gd = ew.global_data

class constants:
    shared_state = {}
    def __init__(self):
        self.__dict__ = self.shared_state

c = constants()
for i in dir(ew):
    if not i.startswith("E_"): continue
    setattr(c, i[2:], getattr(ew, i))

def lgrids(start_grid):
    temp = start_grid
    while temp != None:
        yield start_grid
        temp = start_grid.NextGridThisLevel

def non_block():
    gd.CommunicationDirection = c.COMMUNICATION_POST_RECEIVE
    yield 1
    gd.CommunicationDirection = c.COMMUNICATION_SEND
    yield 2
    ew.CommunicationReceiveHandler()

def EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt):
    Stop = False
    if meta_data.Time >= meta_data.StopTime: Stop = True
    if meta_data.CycleNumber >= meta_data.StopCycle: Stop = True
    meta_data.StartCPUTime = meta_data.CPUTime = LastCPUTime = MPI.Wtime()
    meta_data.LastCycleCPUTime = 0.0

    # Optimized CTP would go here

    # Add Random Forcing if need be
    if gd.RandomForcing:
        for g in lgrids(level_array[0]):
            g.GridData.AppendForcingToBaryonFields()
        Exterior.AppendForcingToBaryonFields()

    # Set top grid boundary condition
    gd.CommunicationReceiveIndex = 0
    gd.CommunicationReceiveCurrentDependsOn = c.COMMUNICATION_NO_DEPENDENCE
    
    for cd in non_block():
        for g in lgrids(level_array[0]):
            g.GridData.SetExternalBoundaryValues(exterior)
            exterior.Prepare(g.GridData)
            ew.CopyOverlappingZones(g.GridData, meta_data, level_array, 0)

    if gd.RandomForcing:
        for g in lgrids(level_array[0]):
            g.GridData.DetachForcingFromBaryonFields()
        Exterior.DetachForcingFromBaryonFields()

top_grid = ew.HierarchyEntry()
meta_data = ew.TopGridData()
exterior = ew.ExternalBoundary()
level_array = ew.LevelHierarchyArray()

ew.CommunicationInitialize([])
#ew.run_enzo_main(["-d", "CollapseTest.enzo"])
gd.LoadBalancing = 0
gd.debug = 1
retval = ew.SetDefaultGlobalValues(meta_data)
restart = False
initial_dt = 0.0

res_fn = "DD0001/moving7_0001"
pf_fn = "CollapseTest.enzo"

if restart:
    ew.Group_ReadAllData(res_fn, top_grid, meta_data, exterior)
    ew.CommunicationPartitionGrid(top_grid, 0)
    #gd.CommunicationDirection = c.COMMUNICATION_SEND_RECEIVE
else:
    initial_dt = ew.InitializeNew(pf_fn, top_grid, meta_data, exterior)

ew.AddLevel(level_array, top_grid, 0)
#ew.EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt)
EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt)
