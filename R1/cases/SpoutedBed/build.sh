
FILE=./mfix-18.1.5  
if [ -e $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   wget https://mfix.netl.doe.gov/source/mfix/mfix-18.1.5.tar.gz
   tar -xzvf mfix-18.1.5.tar.gz 
fi

rsync -av   ../../source/*  $(pwd)  --exclude=usr_rates.f --exclude=cases --exclude=source
MFIX_SRC=$(pwd)/mfix-18.1.5

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles *.d *.mod *.o model post_mfix
   make clean
else
   echo "The directory is good for configuration."
fi


# serial
#debug 
#cmake $MFIX_SRC -DCMAKE_BUILD_TYPE=Debug 
#normal 
echo "configuring mfix-18.1.5 in $MFIX_SRC directory"
cmake $MFIX_SRC -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -Wall -fcheck=all"

#parallel 
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g"
make -j4

echo "run mfixsolver" 
./mfixsolver -f mfix.dat 
