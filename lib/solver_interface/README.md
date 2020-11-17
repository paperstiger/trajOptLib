# How to install

## Dependencies

There are two ways to build the library.
The first way uses CMake and pybind11 (make sure CMake can find
it; do not use pip install pybind11; **install it from source**.)

The second way use setuptools and accepts pybind11 installed from pip (upgrade to the newest version).

Both ways require the installation of at least one of SNOPT and IPOPT.

## Guide

### Obtain SNOPT or IPOPT
However, you might need to set if building SNOPT or IPOPT wrappers.
If you do not build SNOPT wrapper, you have to add `-DBUILD_SNOPT=OFF` to previous command.
Similarly, adding `-DBUILD_IPOPT=OFF` disables IPOPT wrapper.
However, you should install at least one of IPOPT and SNOPT.
To specify SNOPT installation directory, you have to use `-DSNOPT_PATH=path/to/snopt/installation` where include and lib directory exists.
Similarly for IPOPT, you use similar commands `-DIPOPT_PATH=path/to/ipopt/installation`.
You can install IPOPT using `apt-get install coinor-libipopt-dev`, if done so, you could use `-DIPOPT_PATH=/usr` to set IPOPT installation directory.
IPOPT obtained via apt uses a slow linear algebra library and you may consider building it from source to have more control.

### Build using CMake
Use `python setup.py install` to install the package. This requires pybind11 be installed and CMake

In summary, if you installed IPOPT by apt and SNOPT from source at `~/snopt7`, your command is
```bash
python setup.py -DSNOPT_PATH=~/snopt7 -DIPOPT_PATH=/usr install
```
If you install IPOPT by apt and has no SNOPT installed, your command is
```bash
python setup.py -DBUILD_SNOPT=OFF -DIPOPT_PATH=/usr install
```
If you installed SNOPT at `~/snopt7` but has no IPOPT installation, your command is
```bash
python setup.py -DBUILD_IPOPT=OFF -DSNOPT_PATH=~/snopt7 install
```

### Build using setuptools
Use `python setup.py install` to install the package and specify SNOPT, IPOPT directory from commands.

In summary, if you installed IPOPT by apt and SNOPT from source at `~/snopt7`, your command is
```bash
python setup2.py -DSNOPT_PATH=~/snopt7 -DIPOPT_PATH=/usr install
```
If you install IPOPT by apt and has no SNOPT installed, your command is
```bash
python setup2.py -DIPOPT_PATH=/usr install
```
If you installed SNOPT at `~/snopt7` but has no IPOPT installation, your command is
```bash
python setup2.py -DSNOPT_PATH=~/snopt7 install
```


# FAQ
1. What if I passed the wrong arguments to CMake and want to correct it?

CMake use CMakeCache.txt file to store some variables and command line commands sometimes do not change them. In this version, we remove folder `./build` whenever a new installation is performed.

2. What on earth is SNOPT and IPOPT?

Information of SNOPT can be found at <https://web.stanford.edu/group/SOL/snopt.htm>.
IPOPT information is at <https://www.coin-or.org/Ipopt/documentation/>

# Tutorial (coming soon)
