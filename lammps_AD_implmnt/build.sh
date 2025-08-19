#! /bin/bash
mkdir -p build
rm -rf build/*

cmake -B build -S cmake -DPKG_USER-REAXC:BOOL="1" -DWITH_FFMPEG:BOOL="0" -DPKG_MOLECULE:BOOL="1" -DWITH_GZIP:BOOL="0" -DWITH_PNG:BOOL="0" -DBUILD_MPI:BOOL="0" -DWITH_JPEG:BOOL="0"
make --directory=build -j 4
