GSL_HOME=/opt/homebrew/Cellar/gsl/2.7.1
CC = g++
DEBUG = #-g 
GSLINC=-I$(GSL_HOME)/include
GSL=-L$(GSL_HOME)/lib
CFLAGS = -Wall $(DEBUG) $(GSLINC) -O3 -ffast-math -funroll-loops
LFLAGS = -Wall $(DEBUG) $(GSL) -lgsl -lgslcblas -lm


assembly:  tri_tri_intersect.cpp tri_tri_intersect.hpp Geometry.cpp Geometry.hpp MonteCarlo-types.cpp MonteCarlo-types.hpp run_assembly.cpp
	
	$(CC) $(CFLAGS) -c tri_tri_intersect.cpp
	$(CC) $(CFLAGS) -c Geometry.cpp
	$(CC) $(CFLAGS) -c MonteCarlo-types.cpp
	$(CC) $(CFLAGS) -c run_assembly.cpp
	$(CC) $(CFLAGS) tri_tri_intersect.o Geometry.o MonteCarlo-types.o run_assembly.o -o assemble $(LFLAGS)

#fission:  Geometry.cpp Geometry.h MonteCarlo-types.cpp MonteCarlo-types.h run_fission.cpp
	
#	$(CC) $(CFLAGS) -c Geometry.cpp
#	$(CC) $(CFLAGS) -c MonteCarlo-types.cpp
#	$(CC) $(CFLAGS) -c run_fission.cpp
#	$(CC) $(CFLAGS) Geometry.o MonteCarlo-types.o run_fission.o -o fission $(LFLAGS)


clean:
	rm -fr *.o
