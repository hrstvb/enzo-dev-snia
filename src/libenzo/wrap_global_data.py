mappings = {"float":("float", "ebfloat"),
            "double":("double", "double"),
            "int": ("int", "eint")}

output_header = """
# Automatically-generated global data handler
# DO NOT MODIFY DIRECTLY

cdef extern from "../enzo/units.h":
    pass

cdef extern from "../enzo/flowdefs.h":
    pass

cdef extern from "../enzo/CosmologyParameters.h":
    pass

cdef extern from "../enzo/communication.h":
    pass

cdef extern from "../enzo/StarParticleData.h":
    pass

cdef extern from "../enzo/global_data.h":
"""

class_header = """
cdef class _global_data:
"""

var_setter = """
    property %(var)s:
        def __get__(self):
            global %(var)s
            return %(var)s
        def __set__(self, val):
            global %(var)s
            %(var)s = val
"""

loop_setter1 = """
    property %(var)s:
        def __get__(self):
            print "Returning a copy of %(var)s"
            global %(var)s
            retval = []
            for i in range(%(size)s):
                retval.append(%(var)s[i])
            return retval
        def __set__(self, val):
            cdef int i
            global %(var)s
            for i in range(%(size)s):
                %(var)s[i] = val[i]
"""

loop_setter2 = """
    property %(var)s:
        def __get__(self):
            print "Returning a copy of %(var)s"
            global %(var)s
            retval = []
            for i in range(%(size1)s):
                retval.append([])
                for j in range(%(size2)s):
                    retval[-1].append(%(var)s[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            global %(var)s
            for i in range(%(size1)s):
                for j in range(%(size2)s):
                    %(var)s[i][j] = val[i][j]
"""

output_footer = """
# END GLOBAL DATA HANDLER
"""

def parse_global_data():
    defs = " ".join(["-D%s=%s" % (a,b) for a,b in define_macros])
    cmd = "cpp -DPARSE_GLOBAL_DATA %s ../enzo/function_declarations.h > tmp" % defs
    os.system(cmd)
    vars = []
    var_names = []
    for line in (l.strip() for l in open("tmp")):
        if not line.startswith("extern"): continue
        line = line[7:]
        if "//" in line: line = line[:line.rfind("//")]
        if ";" not in line: continue
        if "*" in line: continue
        #if "[" in line: continue
        line = line[:line.rfind(";")]
        c, v = line.split(None, 1)
        if v in var_names: continue
        vars.append((c, v))
        var_names.append(v)

    output = open("enzo_globals.pxi", "w")
    output.write(output_header)
    for ctype, variable in vars:
        if ctype not in mappings: continue
        for var in variable.split(","):
            output.write("    cdef extern %s %s\n" % (
                    ctype.strip(), var.strip()))
    output.write(class_header)
    for ctype, variable in vars:
        if ctype not in mappings: continue
        for var in variable.split(","):
            var = var.strip()
            if "[" in var: 
                if var.count("[") == 1:
                    i1, i2 = var.find("[") + 1, var.find("]")
                    size = var[i1:i2]
                    var = var[:i1-1]
                    output.write(loop_setter1 % dict(var=var,
                                                     size=size))
                if var.count("[") == 2:
                    i1, i2 = var.find("[") + 1, var.find("]")
                    j1, j2 = var.rfind("[") + 1, var.rfind("]")
                    size1, size2 = var[i1:i2], var[j1:j2]
                    var = var[:i1-1]
                    output.write(loop_setter2 % dict(var=var,
                                                     size1=size1,
                                                     size2=size2))
                else: continue
            else:
                output.write(var_setter % dict(var=var))
    output.write(output_footer)

#parse_global_data()