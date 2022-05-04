INSTALL_DIR=$(pwd)
SUNDIALS_DIR=${INSTALL_DIR}/sundials-5.7.0
HYPRE_DIR=${INSTALL_DIR}/hypre-2.24.0
METIS_DIR=${INSTALL_DIR}/metis-5.1.0

#Get Hypre source code
wget -c https://github.com/hypre-space/hypre/archive/refs/tags/v2.24.0.tar.gz -O hypre.tar.gz

#Get Metis 5 source code
wget -c http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz -O metis.tar.gz

#Get SUNDIALS source code
wget -c https://github.com/LLNL/sundials/releases/download/v5.7.0/sundials-5.7.0.tar.gz -O sundials.tar.gz

#Get MFEM source code
git clone https://github.com/mfem/mfem.git

######################
### INSTALL HYPRE
######################
tar -xf hypre.tar.gz
cd ${HYPRE_DIR}/src
./configure --disable-fortran
make -j 2
make install
cd ${INSTALL_DIR}

######################
### INSTALL METIS
######################
tar -xf metis.tar.gz
cd ${METIS_DIR}
make config
make -j 2
mkdir lib
ln -s ../build/Linux-x86_64/libmetis/libmetis.a lib
cd ${INSTALL_DIR}

######################
### INSTALL SUNDAILS
######################
tar -xf sundials.tar.gz
cd ${SUNDIALS_DIR}
mkdir build
cd ${SUNDIALS_DIR}/build
cmake -DCMAKE_INSTALL_PREFIX=${SUNDIALS_DIR}/install -DEXAMPLES_INSTALL_PATH=${SUNDIALS_DIR}/install/examples -DENABLE_MPI:BOOL=ON ${SUNDIALS_DIR}
make install -j 2
cd ${INSTALL_DIR}

######################
### INSTALL MFEM
######################
cd ${INSTALL_DIR}/mfem
mkdir build
cd ${INSTALL_DIR}/mfem/build
cmake -DMFEM_USE_MPI:BOOL=ON -DMFEM_USE_METIS:BOOL=ON -DMFEM_ENABLE_MINIAPPS:BOOL=ON -DMFEM_USE_SUNDIALS:BOOL=ON -DMFEM_USE_ZLIB:BOOL=ON -DHYPRE_DIR=${HYPRE_DIR}/src/hypre -DMETIS_DIR=${METIS_DIR} -DSUNDIALS_DIR=${SUNDIALS_DIR}/install ${INSTALL_DIR}/mfem
make -j 2





