CC=g++
CFLAGS=-O3 -pg

readcif: readcif.cpp unitcell.cpp symmetry.cpp
	$(CC) -o readcif readcif.cpp unitcell.cpp symmetry.cpp $(CFLAGS)

test: testgrid.cpp unitcell.cpp symmetry.cpp energy.cpp
	$(CC) -o testgrid testgrid.cpp energy.cpp unitcell.cpp symmetry.cpp $(CFLAGS)

clean:
	rm readcif
