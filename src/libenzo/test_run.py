import enzo_wrap

#enzo_wrap.run_enzo_main(["-d", "CollapseTest.enzo"])
meta_data = enzo_wrap.TopGridData()
print "MetaData:", meta_data
enzo_wrap.global_data.LoadBalancing = 0
print "LoadBalancing:", enzo_wrap.inquire_LoadBalancing()
retval = enzo_wrap.SetDefaultGlobalValues(meta_data)
print "RetVal:", retval
print "LoadBalancing:", enzo_wrap.inquire_LoadBalancing()
enzo_wrap.global_data.LoadBalancing = 0
print "LoadBalancing:", enzo_wrap.inquire_LoadBalancing()
i = enzo_wrap.global_data.RefineRegionLeftEdge
print i
