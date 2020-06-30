

#wget https://mfix.netl.doe.gov/s3/ef9c88b3/f4a982b34556d34ad702a5dcca5dc337//source/mfix/mfix-20.1.2.tar.gz 
#tar -xzvf mfix-20.1.2.tar.gz 

#if [ -e ./CMakeCache.txt ]; then
#   rm -rf CMakeFiles cmake_install.cmake CMakeFiles
#   make clean
#else
#   echo "The directory is good for mfix compilation."
#fi

echo "remove the files generated in previous run"
rm -rf CASE1.* 

MFIX_SRC=/Users/vkotteda/CODES/mfix-rad/mfix-20.1.2 
#MFIX_SRC=$(pwd)/mfix-20.1.2 

echo "configure mfix"
cmake $MFIX_SRC -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O0 -g -Wall  -Wno-tabs -fcheck=all" -DENABLE_MPI=1
echo "build mfix"
make -j4
echo "run the case"
./mfixsolver -f mfix.dat 
