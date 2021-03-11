## Detailed instructions to install deal.II prerequisities 

### P4EST

1. Download p4est tarball and setup script
```
cd $HOME
mkdir p4est ; cd p4est
wget https://p4est.github.io/release/p4est-2.2.tar.gz
wget https://www.dealii.org/current/external-libs/p4est-setup.sh
``` 

2. Create p4est installation directory, install it and set the appropriate environment variable equal to the installation path

```
cd $HOME
mkdir p4est-install
cd $HOME/p4est
chmod u+x p4est-setup.sh
./p4est-setup.sh p4est-2.2.tar.gz $HOME/p4est_install
export P4EST_DIR=$HOME/p4est_install 
```

### TRILINOS

1. Download trilinos tarball and extract it 

```
cd $HOME
wget https://github.com/trilinos/Trilinos/archive/trilinos-release-13-0-1.tar.gz
tar -xzvf trilinos-release-13-0-1.tar.gz
```  

2. Create installation directory, build directory and install trilinos

```
mkdir Trilinos_install
cd Trilinos-trilinos-release-13-0-1
mkdir build
cd build 
cmake                                                \
    -DTrilinos_ENABLE_Amesos=ON                      \
    -DTrilinos_ENABLE_Epetra=ON                      \
    -DTrilinos_ENABLE_EpetraExt=ON                   \
    -DTrilinos_ENABLE_Ifpack=ON                      \
    -DTrilinos_ENABLE_AztecOO=ON                     \
    -DTrilinos_ENABLE_Sacado=ON                      \
    -DTrilinos_ENABLE_Teuchos=ON                     \
    -DTrilinos_ENABLE_MueLu=ON                       \
    -DTrilinos_ENABLE_ML=ON                          \
    -DTrilinos_ENABLE_ROL=ON                         \
    -DTrilinos_ENABLE_Tpetra=ON                      \
    -DTrilinos_ENABLE_COMPLEX_DOUBLE=ON              \
    -DTrilinos_ENABLE_COMPLEX_FLOAT=ON               \
    -DTrilinos_ENABLE_Zoltan=ON                      \
    -DTrilinos_VERBOSE_CONFIGURE=OFF                 \
    -DTPL_ENABLE_MPI=ON                              \
    -DBUILD_SHARED_LIBS=ON                           \
    -DCMAKE_VERBOSE_MAKEFILE=OFF                     \
    -DCMAKE_BUILD_TYPE=RELEASE                       \
    -DCMAKE_INSTALL_PREFIX=$HOME/Trilinos_install    \
    ../
make install
``` 
Alternatively, the installation process can be completed using multiple processors via

```
make -jN install 
```

where 'N'above is a whole number referring to the number of MPI processes.

