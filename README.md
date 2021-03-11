## Numerical Solution Of Reaction-Diffusion PDEs

Reaction-diffusion(RD) systems form a class of PDEs used to model pattern-formation. In their most basic form, they involve an interplay between
a linear model of diffusion and non-linear model of reaction between multiple species, which can result in a heterogeneous solution(pattern)
under seemingly inconspicuous initial conditions. They arose in the then nascent field of mathematical biology with the foundational work
of Alan Turing, and since then, have been applied to a wide range of scenarios.

This repository is intended to be a collection of C++ codes to model different RD systems. The numerical implementation constitutes a Finite Element(FE) 
approximation to the system of PDEs and the codes are built on the deal.II open source FE library.

## Prerequisities 
FE Package : [**deal.II**](https://www.dealii.org/) **version 9.2**

deal.II prerequisites : [**MPI**](https://www.mpi-forum.org/mpi-40/), [**P4EST**](https://www.p4est.org/), [**BLAS/LAPACK**](https://www.netlib.org/lapack/lug/node11.html), [**Trilinos**](https://trilinos.github.io/)

## deal.II Installation instructions

1. **MPI** : On an HPC platform there will most likely be an MPI implementation which can be loaded as a module. On the desktop the appropriate package manager can be used to install the required libraries. For example, an OpenMPI distribution can be installed in ubuntu as follows :
```
sudo apt-get install libopenmpi-dev openmpi-bin
```

2. **P4EST** : Installation instructions can be found [here](https://www.dealii.org/current/external-libs/p4est.html).

3. **BLAS/LAPACK** : On an HPC platform, there will most likely be some BLAS/LAPACK implementation which can be loaded as a module. On the desktop the appropriate package manager can be used to install the relevant libraries. For example, one can install BLAS/LAPACK in ubuntu as follows :

```
sudo apt-get install libblas3
```
 
4. **Trilinos** : Installation instructions can be found [here](https://www.dealii.org/current/external-libs/trilinos.html).

5. **deal.II** : Once the prerequisites have installed successfully, **deal.II version 9.2.0** can be cloned an installed following the instructions [here](https://www.dealii.org/current/readme.html).

6. After installing **deal.II** make sure to initialize the environment variable *DEAL_II_DIR* as the absolute path to the **deal.II** installation. 

## Compiling the RD code
The code can be compiled as :
```
cmake . ; make 
```
and run as :
```
mpirun -np $NPROCS $EXECUTABLE $INPUTFILE
```

where $NPROCS refers to the number of processors, $EXECUTABLE refers to the problem executable created after compilation, and $INPUTFILE is the input file supplied by the user. For example, in the **Schnakenberg** folder the executable is **RD_PBC** which is generated after compiling **RD_PBC.cc**. The file **testinp.prm** is the input file which is typically modified by the user to provide the required inputs to the problem. 




