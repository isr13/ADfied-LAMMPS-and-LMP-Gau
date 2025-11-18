#!/bin/zsh

for filename in ./*ts*.xyz;
do
        echo ${filename}
        obabel -ixyz ${filename:t:r}.xyz -ogjf -xf format_min.txt -O ${filename:t:r}_min.gjf
done
