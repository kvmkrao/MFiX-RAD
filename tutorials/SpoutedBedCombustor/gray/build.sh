
FILE=./mfix-19.3.0  
if [ -e $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   wget https://mfix.netl.doe.gov/s3/ef9c88b3/1101ec28a35b15cd2d4c94aff939b5be//source/mfix/mfix-19.3.0.tar.gz 
   tar -xzvf mfix-19.3.0.tar.gz ;  
fi

SIZE=$(du -sb mfix-19.3.0 | cut -f1)
   if [[ $SIZE -lt 1000000 ]]; then
      echo "Error in downloading the source file" 
   fi

rsync -av   ../../../source/*  $(pwd)  --exclude=usr_* --exclude=cases 
MFIX_SRC=$(pwd)/mfix-19.3.0

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles *.d *.mod *.o model post_mfix
   make clean
else
   echo "The directory is good for configuration."
fi


echo "configuring mfix-19.3.0 in $MFIX_SRC directory"
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g" 
make -j8

echo "run mfixsolver" 
./mfixsolver -f mfix.dat 
