from distutils.core import setup
from pkgutil import walk_packages


setup(
        name='optsolver',
        version='0.3.0',
        author='Gao Tang',
        author_email='gao.tang@duke.edu',
        #packages=list(find_packages(libpysnopt.__path__, libpysnopt.__name__)),
        scripts=[],
        url='',
        license='LICENSE.txt',
        description='A unified wrapper for many optimization problems',
        #long_description=open('README.md').read(),
        install_requires=[
            'numpy>=1.15.0'
        ],
)
