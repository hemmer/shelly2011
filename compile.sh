#!/bin/bash


# build library
cd ~/code/fourier-methods/
make clean && make
ar -cvq libffs.a fourier.o && echo "Library successfully built!" || (echo "Library failed to build!" && exit -1)

echo

cd ~/code/shelly2011/
make clean && make
