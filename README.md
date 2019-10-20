# trajoptlib
A library for trajectory optimization using direct transcription approach.

##  Installation

Make sure you have

- CMake with version higher than 3.0
- Eigen3 installed somewhere (make sure CMake can find it)
- pybind11 install and cmake can find it (**Do not use pip install**. Use CMake to install from source. For Ubuntu 18.04, you can use `apt-get install python-pybind11` to install it.)

Installation is simple.

0. You have to install pyoptsolver, which is located in `lib/solver_interface`; It is a separate tool for
   defining and solving optimization problems in Python using state-of-the-art nonlinear solvers. Please see instructions there.
1. Install requirements by `pip install -r requirements.txt`
2. Install package using pip, by `pip install -e .` or `python setup.py install`

## Test
You can test a few simple examples. This is done by simply navigating to `test` folder in the main folder. And use `python simpleDemo.py -grad`. You can also explore other possible arguments.
There is a rotor example which requires the following procedure:

1. Navigate to current directory. Make a directory for building `mkdir build & cd build`
2. Use standard CMake install commands: `cmake .. & make & cd ..`

If you did not perform step 1 and 2 in the installation guide, you won't be able to run the rotor example.

## Documentation

Check <https://paperstiger.github.io/trajOptLibDoc/index.html> for documentations.
