#! /bin/bash

# add include to all headers
#FILES=$(ls -1 src/*.h src/USER-REAXC/*.h | grep -v "src/ad_defines.h" | grep -v "src/version.h")
#i=1
#for f in $FILES; do
 # echo "running sed on $f"

  # DEFINE="#ifndef STCE_INC_$i\n  #define STCE_INC_$i\n  #include \"ad_defines.h\"\n#endif\n"
 # DEFINE="#include \"ad_defines.h\"\n"
  #sed -i "1s;^;$DEFINE;" $f
  #i=$((i+1))
#done

# change all double's to ADtype in headers and sources
FILES_C_H=$(ls -1 src/PYTHON/*.h src/PYTHON/*.cpp | grep -v "src/ad_defines.h")
for f in $FILES_C_H; do
  echo "running sed on $f"
  sed -i "s/double/myScalar/g"   $f
  sed -i "s/MY_MPI_Allreduce/MPI_Allreduce/g"   $f
done
