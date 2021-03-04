# Numerical Implementation Of Reaction-Diffusion System

Reaction-diffusion(RD) systems form a class of PDEs used to model pattern-formation. In their most basic form, they involve an interplay between
a linear model of diffusion and non-linear model of reaction between multiple species, which results in a heterogeneous solution(pattern)
under seemingly inconspicuous initial conditions. They arose in the then nascent field of mathematical biology with the foundational work
of Alan Turing, and since then, have been applied to a wide range of scenarios.

This repository is intended to be a collection of C++ codes to model different RD systems, written on top of the deal.II FE library. 

## Prerequisities 
FE Package : [**deal.II version 9.2.0**](https://www.dealii.org/)

deal.II prereqs : **MPI**(Parallel Processing), **P4EST**(Domain partitioning), **BLAS/LAPACK**(Linear Algebra), **Trilinos**(Parallel Linear Algebra)

## deal.II Installation instructions

1. **MPI** : If you are working on some HPC platform, there will most likely be an MPI implementation which can be loaded as a module.
If you are working locally on your desktop, the appropriate package manager can be used to install the relevant libraries. For example, 
in ubuntu, one can install an OpenMPI distribution using 

<div align="center"><i>sudo apt-get install libopenmpi-dev openmpi-bin</i></div>

2. **P4EST** : Installation instructions can be found [here](https://www.dealii.org/current/external-libs/p4est.html)

3. **BLAS/LAPACK** : If you are working on some HPC platform, there will most likely be an BLAS/LAPACK implementation which can be loaded as a module.
If you are working locally on your desktop, the appropriate package manager can be used to install the relevant libraries. For example,
in ubuntu, one can install BLAS/LAPACK using

<div align="center"><i>sudo apt-get install libblas3</i></div>
 
4. **Trilinos** : Installation instructions can be found [here](https://www.dealii.org/current/external-libs/trilinos.html)

## Running the code
The code can be run as follows :

<div align="center"><i>cmake . ; make -j8 ; mpirun -np 4 $EXECUTABLE$ $INPUT FILE$</i></div>

For example, in the **Schnakenberg** folder the executable is **RD_PBC** which is generated after compiling **RD_PBC.cc**. **testinp.prm** is the input file
which is used by the user to provide the required inputs to the problem. 




