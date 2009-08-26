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

for line in open(sys.argv[-1]):
    if not line.startswith(" "*8): continue
    if len(line.strip()) == 0: continue
    if "[" in line:
        skipped.append(line)
        continue
    ls = line.split()
    if ls[0] not in known_types: continue
    for var in ls[1:]:
        print prop_template % dict(var=var, ctype=ls[0])

print "\n" * 2

for skip in skipped:
    print "# SKIPPED: '%s'" % (skip.strip())
