# setup.py

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import os, sys

import numpy
numpy_include = numpy.get_include()

# Let's try to figure out the -D options

for line in open("../enzo/temp.show-flags"):
    if line.startswith("DEFINES"): break

define_macros = []
for opt in line.split("=", 1)[1].split():
    vd = opt[2:].split("=")
    if len(vd) == 1: vd.append(None)
    define_macros.append(tuple(vd))

print define_macros

def check_for_hdf5():
    if "HDF5_DIR" in os.environ:
        return os.environ["HDF5_DIR"]
    elif os.path.exists("hdf5.cfg"):
        return open("hdf5.cfg").read().strip().rstrip()
    print "Reading HDF5 location from hdf5.cfg failed."
    print "Please place the base directory of your HDF5 install in hdf5.cfg and restart."
    print "(ex: \"echo '/usr/local/' > hdf5.cfg\" )"
    sys.exit(1)

H5dir = check_for_hdf5()

setup(name="enzo",
      cmdclass = {'build_ext': build_ext},
      ext_modules=[Extension("enzo_wrap",
                     ["enzo_wrap.pyx"],
                     include_dirs=["/usr/include/", "../enzo/", numpy_include, 
                                   os.path.join(H5dir,"include")],
                     language="c++",
                     libraries=['enzo_p8_b8','hdf5'],
                     library_dirs=["../enzo/",
                                   os.path.join(H5dir,"lib")],
                     define_macros=define_macros
                  )]
      
)
