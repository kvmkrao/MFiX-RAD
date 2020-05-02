
rm -rf CMakeCache.txt  CMakeFiles/ Make* 
MFIX_SRC=/home/vkotteda/CODES/MFiX-RAD/mfix-20.1.2 
#MFIX_SRC=$HOME/Software/mfix_src/mfix-18.1.5
#cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g"
cmake $MFIX_SRC -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O0 -g -Wall  -Wno-tabs -fcheck=all"
make -j4

