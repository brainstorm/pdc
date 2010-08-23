#!/bin/csh
foreach n (`seq 1 1 8`)
    env OMP_NUM_THREADS=$n ./shwater2d
end
