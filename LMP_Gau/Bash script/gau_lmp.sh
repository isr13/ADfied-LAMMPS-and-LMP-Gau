#!/bin/bash

#This script was edited by Indu Sekhar Roy (RWTH Aachen). The original format was for xtb.sh script by Dr. Tian Lu at Beijing Kein Research Center for Natural Sciences (www.keinsci.com)


read atoms derivs charge spin < $2

#Create temporary .xyz file
#the element index should be replaced with element name, and the coordinate should be converted to Angstrom
echo "Generating mol.tmp"
cat >> mol.tmp <<EOF
$atoms

$(sed -n 2,$(($atoms+1))p < $2 | cut -c 1-72)
EOF

echo "Generating mol.xyz via genxyz"
./genxyz

python write_data.py mol.xyz

echo "data file is written, LAMMPS is being called"

rm -f mol.tmp

Location_to_LAMMPS/build/lmp_serial -in in.run


echo "lammps running finished!"

echo "Extracting data from lammps outputs via extderi with derivs= " $derivs
./extderi $3 $atoms $derivs
