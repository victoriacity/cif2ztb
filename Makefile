CC=g++
CFLAGS=-O

readcif: readcif.cpp unitcell.cpp symmetry.cpp
	$(CC) -o readcif readcif.cpp unitcell.cpp symmetry.cpp $(CFLAGS)

clean:
	rm readcif
