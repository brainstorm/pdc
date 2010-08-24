#!/bin/sh

module add summer/10
OUTFILE="blas.out"
PROGRAMFILE="call_blas"

SIZE="400 800 1200 1600 2400 3200 4800 9600"
for k in $SIZE
do
  ./$PROGRAMFILE $k 1>> $OUTFILE
done
