
FILE=./mfix-19.3.0  
if [ -e $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   wget mfix.netl.doe.gov/s3/ef9c88b3/c7a56bea4bf9fc5ae9cfc7a481f4ee06//source/mfix/mfix-19.3.0.tar.gz 
   tar -xzvf mfix-19.3.0.tar.gz ;
fi

SIZE=$(du -sb mfix-19.3.0 | cut -f1)
if [[ $SIZE -lt 1000000 ]]; then
   echo "Error in downloading the source file"
fi

rsync -av   ../../../source/*  $(pwd)  --exclude=usr_rates.f --exclude=cases
MFIX_SRC=$(pwd)/mfix-19.3.0

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles *.d *.mod *.o model post_mfix
   make clean
else
   echo "The directory is good for configuration."
fi


#serial
#debug 
#cmake $MFIX_SRC -DCMAKE_BUILD_TYPE=Debug 
#normal 

echo "configuring mfix-19.1.2 in $MFIX_SRC directory"
#cmake $MFIX_SRC -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -Wall -fcheck=all"
#parallel 
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g"
make -j4

echo "run mfixsolver" 
./mfixsolver -f mfix.dat 
