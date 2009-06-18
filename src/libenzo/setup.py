# setup.py

from distutils.core import setup, Extension
import os

import numpy
numpy_include = numpy.get_include()

setup(name="enzo",
      py_modules=['enzo'], 
      ext_modules=[Extension("_enzo",
                     ["enzo.i"],
                     swig_opts=['-c++'],
                     include_dirs=["../enzo/", numpy_include],
                     libraries=['enzo_p8_b8'],
                     library_dirs=["../enzo/"],
                  )]
      
)
