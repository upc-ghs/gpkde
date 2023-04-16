# gpkde
Fortran code for Grid Projected Kernel Density Estimation of Discrete Particles Distributions

## Overview
The program performs Grid Projected Kernel Density Estimation (GPKDE) of a discrete dataset in one, two or three dimensional domains and is parallelized with the OpenMP library. 

It works as a standalone software by reading an input simulation file, which configures the loading of a source file with data points and additional parameters for defining the reconstruction grid and the optimization for bandwidth selection.

## Command line interface
Some basic command line arguments have been implemented in order to control program execution. These can be requested as help with the instruction ``gpkde --help`` or ``gpkde -h``, which displays the following message in console:
 
```
GPKDE version 0.0.1               
Program compiled Apr 12 2023 19:44:24 with GFORTRAN compiler (ver. *.*.*)       

Fortran code for Grid Projected Kernel Density Estimation of discrete particles distributions

usage:

  gpkde [options] simfile

options:

  -h         --help                Show this message                             
  -l  <str>  --logname    <str>    Write program logs to <str>                   
  -nl        --nolog               Do not write log file                         
  -np <int>  --nprocs     <int>    Run with <int> processes                      
  -p         --parallel            Run in parallel                               
  -v         --version             Show program version                          

For bug reports and updates, follow:                                             
  https://github.com/upc-ghs/gpkde    
```

## Simulation file
For details about the configuration of input files please refer to the program [Documentation](doc/gpkde_IO_001_tmp.pdf). 

## Examples 
A set of possible use cases of the reconstruction module are included as example files in this repository. The currently available simulations include:

- [1D Gaussian distribution](examples/ex01_1dnormal/)
- [2D Heterogeneous distribution](examples/ex02_2dhet/)

## Compilation 
Program includes a `make/Makefile` for the `gfortran` compiler.

## License
MIT License

## Resources
* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)
* [MIT License](https://mit-license.org/)
* [modpath-rw](https://github.com/upc-ghs/modpath-rw)
