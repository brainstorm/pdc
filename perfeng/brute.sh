#!/bin/sh

module add summer/10
OUTFILE="brute1.out"
PROGRAMFILE="scale1"

ARY1="16 32 64 128"
ARY2="128 256 512 1024 2048"
SIZE="400 800 1200 1600 2400 3200 4800 9600"

for j in $ARY2
do
	for i in $ARY1
	do
		echo "BLOCK1 = $i BLOCK2 = $j" 1>> $OUTFILE
		icc -O3 -fast -std=c99 -openmp main.c -o $PROGRAMFILE -DBLOCK1=$i -DBLOCK2=$j
    for k in $SIZE
    do
		  ./$PROGRAMFILE $k 1>> $OUTFILE
    done
		#gcc -O3 -std=c99 main.c -o gscale -lm -DBLOCK1=$i -DBLOCK2=$j
	done
done
