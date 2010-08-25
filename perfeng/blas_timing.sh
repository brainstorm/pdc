#!/bin/sh

module add summer/10
OUTFILE="brute_temp.out"
PROGRAMFILE="scale"

SIZE="400 800 1200 1600 2400 3200 4800 9600"
for k in $SIZE
do
  ./$PROGRAMFILE $k 1>> $OUTFILE
done
