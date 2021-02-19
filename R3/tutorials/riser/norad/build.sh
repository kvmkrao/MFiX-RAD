
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

MFIX_SRC=$(pwd)/mfix-19.3.0

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles *.d *.mod *.o model post_mfix
   make clean
else
   echo "The directory is good for configuration."
fi


echo "configuring mfix-19.3.0 in $MFIX_SRC directory"
#serial
#cmake $MFIX_SRC -DENABLE_MPI=0 -DCMAKE_Fortran_FLAGS="-O2 -g"
#parallel 
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g" 
make -j8

echo "run mfixsolver" 
./mfixsolver -f mfix.dat 
