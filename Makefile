CC=g++
NVCC=nvcc
CFLAGS=-O3 -pg
CUDAFLAGS=-O3 --use_fast_math
readcif: readcif.cpp unitcell.cpp symmetry.cpp
	$(CC) -o readcif readcif.cpp unitcell.cpp symmetry.cpp $(CFLAGS)

test: testgrid.cpp unitcell.cpp symmetry.cpp energy.cpp
	$(CC) -o testgrid testgrid.cpp energy.cpp unitcell.cpp symmetry.cpp $(CFLAGS)

testcuda: testgrid.cu unitcell.cpp symmetry.cpp energy.cu
	$(NVCC) -o testgrid_cuda testgrid.cu energy.cu unitcell.cpp symmetry.cpp $(CUDAFLAGS)

clean:
	rm readcif
