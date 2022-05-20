cd ~/git/FMM3D/
make install
cd ~/git/01_TES_Package/Apr20_2022_fmm3dbie
gfortran -O3 -fPIC -fopenmp -march=native -std=legacy convert_geom_go3.f90 -lfmm3d -lfmm3dbie -L/mnt/home/mrachh/lib
./a.out
