import os
import re
import sys
from glob import glob
import sysconfig
import platform
import subprocess
import shutil

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from pybind11.setup_helpers import Pybind11Extension


files = ["wrapper.cpp"]
SNOPTWRAPPERSRC = ["snopt-interface/src/funcStyle.cpp", "snopt-interface/src/toyfunction.cpp"]
IPOPTWRAPPERSRC = ["ipopt-interface/src/ipoptWrapper.cpp"]

inc_dirs = ["common", ""]
lib_dirs = []
link_libs = []
defs = []
has_snopt = False
has_ipopt = False
newargs = []
for arg in sys.argv:
    if arg.startswith('-D'):
        choice, direc = arg.split('=')
        inc_dirs.append(os.path.join(direc, "include"))
        lib_dirs.append(os.path.join(direc, "lib"))
        if 'SNOPT' in choice:
            link_libs.extend(['snopt7', 'snopt7_cpp'])
            inc_dirs.append("snopt-interface/include")
            files.extend(SNOPTWRAPPERSRC)
            defs.append(("SNOPT", None))
            has_snopt = True
        if 'IPOPT' in choice:
            link_libs.append('ipopt')
            inc_dirs.append("ipopt-interface/include")
            files.extend(IPOPTWRAPPERSRC)
            defs.append(("IPOPT", None))
            defs.append(("ENABLEIP", None))
            defs.append(("HAVE_CSTDDEF", None))
            has_ipopt = True
    else:
        newargs.append(arg)
sys.argv = newargs

if not (has_snopt or has_ipopt):
    raise Exception("At least one of SNOPT and IPOPT has to be installed and specified by -DSNOPT= or -DIPOPT=")

print(files)
print(inc_dirs)
print(link_libs)
print(lib_dirs)

# add stuff to them
base = "pyoptsolver/src"
files = [os.path.join(base, file_) for file_ in files]
inc_dirs = [os.path.join(base, inc_) for inc_ in inc_dirs]
lib_dirs = [os.path.join(base, lib_) for lib_ in lib_dirs]
if os.path.exists('./build'):
    shutil.rmtree('./build')

ext_modules = [
        Pybind11Extension(
            "pyoptsolvercpp",
            files,
            include_dirs=inc_dirs,
            libraries=link_libs,
            library_dirs=lib_dirs,
	    define_macros=defs
            )
        ]

setup(
    name='pyoptsolver',
    version='0.7.0',
    author='Gao Tang',
    author_email='gaotang2@illinois.edu',
    license='LICENSE.txt',
    description='A unified wrapper for many optimization problems',
    ext_modules=ext_modules,
    packages=find_packages(),
    package_dir={'': './'},
    zip_safe=False,
)
