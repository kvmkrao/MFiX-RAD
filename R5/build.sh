

wget https://mfix.netl.doe.gov/s3/ef9c88b3/19907d38a85d175eee92f114b1b4d085//source/mfix/mfix-20.4.2.tar.gz
tar -xzvf mfix-20.4.2.tar.gz 

if [ -e ./CMakeCache.txt ]; then
   rm -rf CMakeFiles cmake_install.cmake CMakeFiles
   make clean
else
   echo "The directory is good for mfix compilation."
fi

echo "remove the files generated in previous run"
rm -rf CASE1.* 

MFIX_SRC=$(pwd)/mfix-20.4.2 

#MFIX_SRC=/home/vkotteda/Software/mfix_src/mfix-20.2.0/
MFIX_SRC=/Users/vkotteda/Software/mfix-rad/mfix-20.4.2

echo "configure mfix"
#cmake $MFIX_SRC -DCMAKE_Fortran_FLAGS="-O0 -g -Wall  -Wno-tabs -fcheck=all -fbacktrace"
#-ffpe-trap=invalid,zero,overflow,underflow,denormal"  
cmake $MFIX_SRC -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O0 -g -Wall  -Wno-tabs -fcheck=all" -DENABLE_MPI=1
echo "build mfix"
make -j4
#echo "run the case"
#./mfixsolver -f mfix.dat 
