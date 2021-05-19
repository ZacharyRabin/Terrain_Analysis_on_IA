from setuptools import setup

from Cython.Build import cythonize

from distutils.core import setup, Extension

import glob

cpp_files = glob.glob('**/*.cpp', recursive=True)

cpp_files.remove('main.cpp')

if 'PyMain.cpp' in cpp_files:
    cpp_files.remove('PyMain.cpp')


ext = [Extension("PythonMain",
    sources=["PyMain.pyx"] + cpp_files,
    language="c++",
    include_dirs=['..','../core_library/sources','core_library/sources/'])]
setup(ext_modules=cythonize(ext))

# use command 
# python setup.py build_ext --inplace
# to compile and create wrapper