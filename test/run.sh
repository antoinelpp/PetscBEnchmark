#!/usr/bin/bash

cp ../build/src/lppic .
time mpirun -n 2 ./lppic -ksp_type cg  -pc_type gamg -xmax 256 -ymax 1000 -log_view 1
time mpirun -n 4 ./lppic -ksp_type bcgs  -pc_type hypre -xmax 256 -ymax 2000
