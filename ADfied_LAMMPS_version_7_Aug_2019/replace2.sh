#! /bin/bash

# change all double's to ADtype in headers and sources
FILES_C_H=$(ls -1 src/*.h src/*.cpp src/USER-REAXC/*.h src/USER-REAXC/*.cpp | grep -v "src/ad_defines.h")
for f in $FILES_C_H; do
  echo "running sed on $f"
  sed -i "s/MY_MPI_Allreduce/MPI_Allreduce/g"   $f

done
