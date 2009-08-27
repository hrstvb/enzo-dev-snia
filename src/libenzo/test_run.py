from yt.mods import *
import enzo_wrap.enzo_routines as er

vals = er.main(True, 'DD0000/moving7_0000')
pf = load("DD0000/moving7_0000")
top_grid, meta_data, exterior, level_array, initial_dt = vals

g1=level_array[0].GridData
g2=level_array[1].GridData
g3=level_array[1].NextGridThisLevel.GridData
field = g1.get_baryon_field(0)
pvel = g1.ParticleVelocity
pnum = 0
for g in [g1, g2, g3]: pnum += g.ParticleNumber.size
print pf.h.gridNumberOfParticles.sum()
print g1.NumberOfParticles + g2.NumberOfParticles + g3.NumberOfParticles
print pnum
