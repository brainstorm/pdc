#!/bin/sh

ARY="16 32 64 128 256 512 1024 2048"

for j in $ARY
do
	for i in $ARY
	do
		echo "BLOCK1= $i, BLOCK2= $j"
		icc -O3 -fast -std=c99 -openmp main.c -o scale
		#gcc -O3 -std=c99 main.c -o gscale -lm -DBLOCK1=$i -DBLOCK2=$j
		./scale 1000
	done
done