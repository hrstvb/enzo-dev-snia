import enzo_wrap.enzo_routines as er

vals = er.main(True, 'DD0000/moving7_0000')
top_grid, meta_data, exterior, level_array, initial_dt = vals

g=level_array[0].GridData
g.get_baryon_field(0)
