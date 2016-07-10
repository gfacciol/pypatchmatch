from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

files = ['patchmatch.pyx', 'img.cc']

extra_compile_args=["-DDONT_USE_MAIN"]
extra_link_args=[]
## if the compiler supports OpenMP uncomment:
#extra_compile_args.append('-fopenmp')
#extra_link_args.append('-fopenmp')

extensions = [Extension("patchmatch", files, 
   language="c++", 
   include_dirs=[np.get_include()], 
   extra_compile_args=extra_compile_args, 
   extra_link_args=extra_link_args)]


setup( ext_modules = cythonize(extensions))
