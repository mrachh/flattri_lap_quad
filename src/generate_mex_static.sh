FMMLIBNAME=-lfmm3d
FMM3DBIELIBNAME=-lfmm3dbie_matlab
FMMPATH=/mnt/home/mrachh/lib
FMM3DBIEPATH=/mnt/home/mrachh/lib
LLNAME=flattrirouts
MWNAME=${LLNAME}.mw
GATEWAY=${LLNAME}.c
rm -rf lpeval_flattri.o
rm -rf potint2.o
gfortran -fPIC -march=native -O3 -c lpeval_flattri.f90 -o lpeval_flattri.o -L${FMMPATH} ${FMMLIBNAME} -L${FMM3DBIEPATH} ${FMM3DBIELIBNAME}
gfortran -fPIC -march=native -O3 -c potint2.f90 -o potint2.o -L${FMMPATH} ${FMMLIBNAME} -L${FMM3DBIEPATH} ${FMM3DBIELIBNAME}
mex -v ${GATEWAY} lpeval_flattri.o potint2.o /mnt/home/mrachh/lib/libfmm3d.a /mnt/home/mrachh/lib/libfmm3dbie_matlab.a -largeArrayDims -DMWF77_UNDERSCORE1 -output ${LLNAME} -L${FMMPATH} -ldl -lm -lstdc++ -lgfortran ${FMMLIBNAME} -lgomp -lblas -llapack
