#! /bin/bash

#wget https://github.com/lammps/lammps/archive/refs/tags/stable_7Aug2019.tar.gz
tar -xzf stable_7Aug2019.tar.gz
#git clone https://github.com/lammps/lammps.git && cd lammps && git checkout stable_7Aug2019

LAMMPS_DIR=lammps-stable_7Aug2019

cp ad.hpp ad_defines.h $LAMMPS_DIR/src
cp build.sh $LAMMPS_DIR

cd $LAMMPS_DIR

# add include if ad_defines.h to all headers which don't already have it
FILES_H=$(ls -1 src/*.h src/USER-REAXC/*.h | grep -v "ad_defines.h" | grep -v "src/version.h")

for f in $FILES_H; do
   echo "running sed on $f"

   DEFINE="#include \"ad_defines.h\"\n"
   sed -i "1s;^;$DEFINE;" $f
done

# change all double's to ADtype in headers and sources
FILES_C_H=$(ls -1 src/*.h src/*.cpp src/USER-REAXC/*.h src/USER-REAXC/*.cpp | grep -v "ad_defines.h")
for f in $FILES_C_H; do
  echo "running sed on $f"
  sed -i "s/double/myScalar/g"   $f
done

# change all calls to MPI_Allreduce with wrapped version MY_MPI_Allreduce
#FILES_C_H=$(ls -1 src/*.h src/*.cpp src/USER-REAXC/*.h src/USER-REAXC/*.cpp | grep -v "src/ad_defines.h")
#for f in $FILES_C_H; do
#  echo "running sed on $f"
#  sed -i "s/MY_MPI_Allreduce/MPI_Allreduce/g" $f
#done