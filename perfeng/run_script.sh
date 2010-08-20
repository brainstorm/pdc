#!/bin/sh 
 for i in `seq 10 1 200`;
do
./scale $i >> output
done
