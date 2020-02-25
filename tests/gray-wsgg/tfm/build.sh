
FILE=./mfix-19.2.1  
if [ -e $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   wget https://mfix.netl.doe.gov/s3/ef9c88b3/d2b4a9d0820a7a0b6099f703e6a295ef//source/mfix/mfix-19.2.1.tar.gz 
   tar -xzvf mfix-19.2.1.tar.gz ;  
fi

rsync -av   ../../../source/*  $(pwd)  --exclude=usr_rates.f --exclude=cases 
MFIX_SRC=$(pwd)/mfix-19.2.1

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles *.d *.mod *.o model post_mfix
   make clean
else
   echo "The directory is good for configuration."
fi


echo "configuring mfix-19.1.2 in $MFIX_SRC directory"
#serial
#cmake $MFIX_SRC -DENABLE_MPI=0 -DCMAKE_Fortran_FLAGS="-O2 -g"
#parallel 
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g" 
make -j8

echo "run mfixsolver" 
./mfixsolver -f mfix.dat 
