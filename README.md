# trajOptLib
A library for trajectory optimization using direct transcription approach.

##  Installation

Make sure you have

- CMake with version higher than 3.0
- Eigen3 installed somewhere
- pybind11 install and cmake can find it

Installation is simple. 

1. In terminal navigate to current directory
2. Create new directory called Build by "mkdir Build" and go inside "cd Build"
3. Invoke CMake to build the library by "cmake ..", followed by "make"
4. Install package using pip, by "pip install -e .."

You can test a few simple examples. This is done by simply navigating to "test" folder in the main folder. And use "python simpleDemo.py -grad". You can also explore other possible arguments.

Check <https://paperstiger.github.io/trajOptLibDoc/index.html> for documentations.

Optionally, you can use Ipopt as the solver for solving problems. 

Check <https://github.com/xuy/pyipopt> for Python binding of the solver <https://projects.coin-or.org/Ipopt>
