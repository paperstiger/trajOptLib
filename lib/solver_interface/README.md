# Installation

Use python setup.py develop to install the package.
However, you might need to set if building SNOPT or IPOPT wrappers.
If you do not build SNOPT wrapper, you have to add -DBUILD_SNOPT=OFF to previous command.
Similarly, add -DBUILD_IPOPT=OFF disables IPOPT wrapper.
To specify SNOPT installation directory, you have to use -DSNOPT_PATH=path/to/snopt/installation where include and lib directory exists.
Similarly for IPOPT, you use similar commands.
You can install IPOPT using apt-get install coinor-libipopt-dev
