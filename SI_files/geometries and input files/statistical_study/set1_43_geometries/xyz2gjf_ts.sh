#!/bin/zsh

for filename in ./*ts*.xyz;
do
        echo ${filename}
        obabel -ixyz ${filename:t:r}.xyz -ogjf -xf format.txt -O ${filename:t:r}.gjf
done
