#!/usr/bin/bash

cd ../build
cmake ..
make
cd ../test

cp ../build/src/lppic .
time mpirun -n 2 ./lppic -ksp_type cg  -pc_type gamg -xmax 256 -ymax 1000 -log_view 0
time mpirun -n 4 ./lppic -lib hypre -xmax 256 -ymax 2000
