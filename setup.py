from distutils.core import setup
from pkgutil import walk_packages

import trajOptLib

def find_packages(path, prefix=""):
    yield prefix
    prefix = prefix + "."
    for _, name, ispkg in walk_packages(path, prefix):
        if ispkg:
            yield name

setup(
        name='trajOptLib',
        version='0.3.0',
        author='Gao Tang',
        author_email='gao.tang@duke.edu',
        packages=list(find_packages(trajOptLib.__path__, trajOptLib.__name__)),
        scripts=[],
        url='',
        license='LICENSE.txt',
        description='A serious tools for trajectory optimization',
        long_description=open('README.md').read(),
        install_requires=[
            'numpy>=1.13.0',
            'autograd',
            'scipy>=1.0.0',
            'matplotlib>=1.13.0'
        ],
)
