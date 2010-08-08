import types, cPickle, h5py, copy
import numpy as na

funcs = (types.BuiltinFunctionType, types.MethodType, types.FunctionType,
         types.BuiltinMethodType)

def serialize(this_grid):
    grid_dict = {}
    for i in dir(this_grid):
        if i.startswith("_"): continue
        f = getattr(this_grid, i)
        if isinstance(f, funcs):
            continue # We will handle these later
        # This should copy arrays and lists of arrays correctly
        grid_dict[i] = copy.deepcopy(f) 
    # Copy baryon fields
    fdata = [[], []]
    fluxes = []
    for ind in range(this_grid.NumberOfBaryonFields):
        # Old then new
        fdata[0].append(this_grid.get_baryon_field(ind, True)[:])
        fdata[1].append(this_grid.get_baryon_field(ind, False)[:])
        fluxes.append(copy.deepcopy(this_grid.get_fluxes(ind)))
    grid_dict["fields"] = fdata
    grid_dict["fluxes"] = fluxes
    return grid_dict

do_first = ["GridRank", "GridDimension"]

def deserialize(grid_dict, this_grid):
    this_grid.GridRank = grid_dict.pop("GridRank")
    this_grid.GridDimension = grid_dict.pop("GridDimension")
    for i in dir(this_grid):
        if i.startswith("_"): continue
        f = getattr(this_grid, i)
        if isinstance(f, funcs):
            continue # We will handle these later
        if isinstance(f, na.ndarray):
            f[:] = grid_dict[i]
        elif f == None:
            if grid_dict[i] != None: print "Skipping ", i
        else:
            setattr(this_grid, i, grid_dict[i])
    # Now 
