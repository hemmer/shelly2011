#!/bin/bash

N=$1; dt=$2; tmax=$3; Wi=$4
outputFolder=`printf "output/N_%03d_Wi%05.02f_dt%g" $N $Wi $dt`
mkdir -p $outputFolder

cp scripts/process_images.sh $outputFolder
cp Makefile $outputFolder
cp shelly.{c,h} $outputFolder

cd $outputFolder

rm -f st-*
rm -f images/*
rm -f plotting_*

make clean && make

echo "./shelly $@ > output-data"
./shelly $@ > output-data

./process_images.sh trS
./process_images.sh vort
./process_images.sh v

./process_images.sh Kxx
./process_images.sh Kyx
./process_images.sh Kxy
./process_images.sh Kyy

exit 0
