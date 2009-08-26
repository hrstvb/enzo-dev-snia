from mpi4py import MPI
import enzo_wrap as ew

class constants:
    shared_state = {}
    def __init__(self):
        self.__dict__ = self.shared_state

c = constants()
for i in dir(ew):
    if not i.startswith("E_"): continue
    setattr(c, i[2:], getattr(ew, i))

top_grid = ew.HierarchyEntry()
meta_data = ew.TopGridData()
exterior = ew.ExternalBoundary()
level_array = ew.LevelHierarchyArray()

ew.CommunicationInitialize([])
#ew.run_enzo_main(["-d", "CollapseTest.enzo"])
ew.global_data.LoadBalancing = 0
ew.global_data.debug = 1
retval = ew.SetDefaultGlobalValues(meta_data)
ew.read_all_data("DD0000/moving7_0000", top_grid, meta_data, exterior)
ew.global_data.CommunicationDirection = c.COMMUNICATION_SEND_RECEIVE
ew.CommunicationPartitionGrid(top_grid, 0)
ew.AddLevel(level_array, top_grid, 0)

ew.EvolveHierarchy(top_grid, meta_data, exterior, level_array, 0.0)
