from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

files = ['patchmatch.pyx']

extensions = [Extension("patchmatch", files, 
   language="c++", include_dirs=[np.get_include()])]

setup( ext_modules = cythonize(extensions))
