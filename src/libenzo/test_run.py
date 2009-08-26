from mpi4py import MPI
import enzo_wrap

#enzo_wrap.run_enzo_main(["-d", "CollapseTest.enzo"])
meta_data = enzo_wrap.TopGridData()
print "MetaData:", meta_data
enzo_wrap.global_data.LoadBalancing = 0
retval = enzo_wrap.SetDefaultGlobalValues(meta_data)
i = enzo_wrap.global_data.RefineRegionLeftEdge
print i
enzo_wrap.read_all_data("DD0000/moving7_0000", meta_data)
