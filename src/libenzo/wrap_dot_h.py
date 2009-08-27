# This accepts a .pxi file as an argument
# For every line in the file that starts with eight spaces and one of the
# types listed below, it generates a proxy property and outputs that.

import sys

known_types = set(["Eflt", "FLOAT", "Eint", "float", "double", "char *"])

skipped = []

prop_template = """
    property %(var)s:
        def __get__(self):
            return self.thisptr.%(var)s
        def __set__(self, %(ctype)s val):
            self.thisptr.%(var)s = val
"""

loop_setter1 = """
    property %(var)s:
        def __get__(self):
            print "Returning a copy of %(var)s"
            retval = []
            for i in range(%(size)s):
                retval.append(self.thisptr.%(var)s[i])
            return retval
        def __set__(self, val):
            cdef int i
            for i in range(%(size)s):
                self.thisptr.%(var)s[i] = val[i]
"""

loop_setter2 = """
    property %(var)s:
        def __get__(self):
            print "Returning a copy of %(var)s"
            retval = []
            for i in range(%(size1)s):
                retval.append([])
                for j in range(%(size2)s):
                    retval[-1].append(self.thisptr.%(var)s[i][j])
            return retval
        def __set__(self, val):
            cdef int i, j
            for i in range(%(size1)s):
                for j in range(%(size2)s):
                    self.thisptr.%(var)s[i][j] = val[i][j]
"""


for line in open(sys.argv[-1]):
    if not line.startswith(" "*8): continue
    if len(line.strip()) == 0: continue
    ls = line.split()
    if ls[0] not in known_types: continue
    for var in ls[1:]:
        if "*" in var:
            skipped.append((1, ls[0], var, line))
            continue
        elif "(" in var or ")" in var: continue
        var = var.strip()
        if "[" in var:
            if var.count("[") == 1:
                i1, i2 = var.find("[") + 1, var.find("]")
                size = var[i1:i2]
                var = var[:i1-1]
                print loop_setter1 % dict(var=var,
                                          size=size)
            elif var.count("[") == 2:
                i1, i2 = var.find("[") + 1, var.find("]")
                j1, j2 = var.rfind("[") + 1, var.rfind("]")
                size1, size2 = var[i1:i2], var[j1:j2]
                var = var[:i1-1]
                print loop_setter2 % dict(var=var,
                                          size1=size1,
                                          size2=size2)
            else: skipped.append((2, ls[0], var.strip(), line))
        else:
            print prop_template % dict(var=var, ctype=ls[0])

print "\n" * 2

for skip in skipped:
    print "# SKIPPED: '%s'" % (skip[2])
