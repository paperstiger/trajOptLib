# trajOptLib
A library for trajectory optimization using direct transcription approach.

##  Installation

Make sure you have

- CMake with version higher than 3.0
- Eigen3 installed somewhere
- pybind11 install and cmake can find it

Installation is simple. 

0. You have to install pyoptsolver, which is located in lib/solver_interface; see instructions there.
1. Navigate to current directory. Make a directory for building "mkdir build & cd build"
2. Use standard CMake install commands: "cmake .. & make & cd .."
3. Install requirements by "pip install -r requirements.txt"
3. Install package using pip, by "pip install -e .." or "python setup.py install"

You can test a few simple examples. This is done by simply navigating to "test" folder in the main folder. And use "python simpleDemo.py -grad". You can also explore other possible arguments.

Check <https://paperstiger.github.io/trajOptLibDoc/index.html> for documentations.

Optionally, you can use Ipopt as the solver for solving problems. 

Check <https://github.com/xuy/pyipopt> for Python binding of the solver <https://projects.coin-or.org/Ipopt>
