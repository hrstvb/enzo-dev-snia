# setup.py

#from distutils.core import setup, Extension
#from Cython.Distutils import build_ext
from mpidistutils import setup
from mpidistutils import Distribution, Extension, Executable
from mpidistutils import config, build, install, clean
from mpidistutils import build_ext, build_exe
from mpidistutils import install_data, install_exe

import os, sys, subprocess

import numpy
numpy_include = numpy.get_include()

# Let's try to figure out the -D options

for line in open("../enzo/temp.show-flags"):
    if line.startswith("DEFINES"): break

define_macros = [("SHARED_LIBRARY",None)]
for opt in line.split("=", 1)[1].split():
    vd = opt[2:].split("=")
    if len(vd) == 1: vd.append(None)
    define_macros.append(tuple(vd))

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

for i,j in define_macros:
    if i.startswith("CONFIG_BFLOAT_"):
        bfloat = i.split("_")[-1]
    if i.startswith("CONFIG_PFLOAT_"):
        pfloat = i.split("_")[-1]

args = ["cython", "-a", "--cplus", "enzo_module.pyx"]
print args
p = subprocess.Popen(args, cwd=os.getcwd() + "/enzo_wrap/")
p.communicate()
if p.returncode:
    print p.returncode
    sys.exit(1)

print "Linking against %s-%s" % (pfloat, bfloat)
setup(name="enzo", packages = "enzo_wrap",
      cmdclass = {'build_ext': build_ext},
      ext_modules=[Extension("enzo_wrap/enzo_module",
                     ["enzo_wrap/enzo_module.cpp"],
                     include_dirs=["/usr/include/", "../enzo/", numpy_include, 
                                   os.path.join(os.getcwd(), "enzo_wrap/"), 
                                   os.path.join(H5dir,"include")],
                     language="c++",
                     libraries=['enzo_p%s_b%s' % (pfloat, bfloat),'hdf5'],
                     library_dirs=["../enzo/",
                                   os.path.join(H5dir,"lib")],
                     define_macros=define_macros
                  )]
      
)
