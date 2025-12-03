#! /bin/bash

# change all calls to MPI_Allreduce with wrapped version MY_MPI_Allreduce
FILES_C_H=$(ls -1 src/*.h src/*.cpp src/USER-REAXC/*.h src/USER-REAXC/*.cpp | grep -v "src/ad_defines.h")
for f in $FILES_C_H; do
  echo "running sed on $f"
  sed -i "s/MPI_Allreduce/MY_MPI_Allreduce/g" $f
done
