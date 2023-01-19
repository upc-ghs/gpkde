# gpkde
Fortran module for Grid Projected Kernel Density Estimation on a regular structured grid. 

## Overview
The module is integrated into `modpath-rw` and is responsible for providing smoothed concentration reconstruction. 
It can also be used as a standalone, and once a list of points `(x,y,z)` is provided, program will output the list of grid nodes with their corresponding histogram and smoothed reconstruction. 

## Compilation 
Program includes a `make/Makefile` for `gfortran` compilers.

### Note
For compilation of `modpath-rw`, this repository needs to be included in the root directory of the random walk model. Makefiles for that project assume this directory is present.


## Resources
* [modpath-rw](https://github.com/upc-ghs/modpath-rw)
* [modpath-omp](https://github.com/upc-ghs/modpath-omp)
* [flopy](https://github.com/modflowpy/flopy)
* [flopyrw](https://github.com/upc-ghs/flopyrw)
