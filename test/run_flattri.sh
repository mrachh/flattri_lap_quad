gfortran -O3 -fPIC -fopenmp -march=native -std=legacy flattriquad.f90 ../src/lpeval_flattri.f90 ../src/potint2.f90 -lfmm3d -lfmm3dbie -L/mnt/home/mrachh/lib
./a.out
