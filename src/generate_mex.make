MEX=mex
FMM3DBIELIBPATH=C:/git/fmm3dbie/lib-static/libfmm3dbie_matlab.a
FMM3DLIBPATH=C:/git/FMM3D-sep20-2021/lib-static/libfmm3d.a
MFLAGS = -largeArrayDims -DMWF77_UNDERSCORE0 -D_OPENMP
MLIBS = -lmex -lmat -lmx -lgfortran -lquadmath -lgomp -lmwblas -lmwlapack -LC:/lib -L"C:/MATLAB/R2020b/bin/win64"
MINGW_PATH=C:/mingw-w64/mingw64/lib/gcc/x86_64-w64-mingw32/6.3.0
LLNAME=flattrirouts
GATEWAY=${LLNAME}.c
BLASEQ=-Ddgemm_=dgemm

FC=gfortran
FFLAGS=-fPIC -O3 -fopenmp -fno-underscoring

.PHONY: all clean

%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

OBJECTS=lpeval_flattri.o potint2.o

all: matlab

matlab: $(OBJECTS)
	$(MEX) -v ${GATEWAY} ${OBJECTS} ${FMM3DBIELIBPATH} ${FMM3DLIBPATH} ${MFLAGS} -output ${LLNAME} ${MLIBS} -L${MINGW_PATH}

clean:
	rm -f $(OBJECTS)
	rm -f *.mex*
