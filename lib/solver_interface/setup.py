import os
import re
import sys
import glob
import sysconfig
import platform
import subprocess
import shutil

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]
        cmake_args.extend(add_args)
        print(cmake_args)

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j8']

        # cmake_args += ["-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"]
        # cmake_args += ["-DCMAKE_INSTALL_RPATH={}".format("$ORIGIN")]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output

add_args = []
for arg in sys.argv:
    if arg.startswith('-D'):
        add_args.append(arg)
sys.argv = list(filter(lambda x: not x.startswith('-D'), sys.argv))
# remove build directory if it exists
if os.path.exists('./build'):
    shutil.rmtree('./build')
print('Arguments passed to CMake:', add_args)
major_ver, minor_ver = sys.version_info[:2]

setup(
    name='pyoptsolver',
    version='0.6.0',
    author='Gao Tang',
    author_email='gaotang2@illinois.edu',
    license='LICENSE.txt',
    description='A unified wrapper for many optimization problems',
    ext_modules=[CMakeExtension('pyoptsolver/pyoptsolvercpp')],
    packages=find_packages(),
    package_dir={'': './'},
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[
        'numpy'
    ],
)
