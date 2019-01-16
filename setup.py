import os
import re
import sys
import sysconfig
import platform
import subprocess
import copy

from pkgutil import walk_packages
from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import trajOptLib


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
            build_args += ['--', '-j2']

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


def find_packages(path, prefix=""):
    yield prefix
    prefix = prefix + "."
    for _, name, ispkg in walk_packages(path, prefix):
        if ispkg:
            yield name


print(sys.argv)
copied = copy.copy(sys.argv)
add_args = []
for arg in sys.argv:
    print(arg)
    if arg.startswith('-D'):
        add_args.append(arg)
        copied.remove(arg)
print(sys.argv)
sys.argv = copied


setup(
        name='trajOptLib',
        version='0.4.0',
        author='Gao Tang',
        author_email='gao.tang@duke.edu',
        packages=list(find_packages(trajOptLib.__path__, trajOptLib.__name__)),
        scripts=[],
        url='',
        license='LICENSE.txt',
        description='A serious tools for trajectory optimization',
        long_description=open('README.md').read(),
        cmdclass=dict(build_ext=CMakeBuild),
        zip_safe=False,
        install_requires=[
            'numpy>=1.15.0',
            'autograd',
            'scipy>=1.0.0',
            'matplotlib>=1.13.0',
            'pyoptsolver'
        ],
)
