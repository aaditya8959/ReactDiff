## Detailed instructions to install deal.II prerequisities 

### P4EST

1. Download p4est tarball and setup script
```
cd $HOME
mkdir p4est ; cd p4est
wget https://p4est.github.io/release/p4est-2.2.tar.gz
wget https://www.dealii.org/current/external-libs/p4est-setup.sh
``` 

2. Create p4est installation directory, install it and set environmental variable to installation path

```
cd $HOME
mkdir p4est-install
cd $HOME/p4est
chmod u+x p4est-setup.sh
./p4est-setup.sh p4est-2.2.tar.gz $HOME/p4est-install
export P4EST_DIR=$HOME/p4est-install 
```



