PyEDA is a Python3 package design to easily perform Essential Dynamics analysis
of proteins. It can be imported to write your own analysis or use the
predefined script (with a GUI or CLI, see below) that calculates the
eigenvectors of the covariance matrix, plots them and projects the structure
onto a certain number of principal eigenvectors and generates a trajectory
corresponding to the eigenvector.


Installation instructions:

To install PyEDA in your computer you just need a Python3 distribution
installed, if you have it just type:

    python3 setup.py install

from the directory where the code is. This will install the package in your
python path and you will be able to import it.

If you just want to run the pre-made script you can use:

    python3 -m PyEDA

to run the graphical interface, or

    python3 -m PyEDA -ng

to use it via the command line, for a complete description of all the options,
type:

    python3 -m PyEDA -h
