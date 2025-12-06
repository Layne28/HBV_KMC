CC = g++
DEBUG = #-g 
#GSLINC=-I$(GSL_HOME)/include
GSLINC=-I$(HOME)/.local/include/gsl
#GSL=-L$(GSL_HOME)/lib
GSL=-L$(HOME)/.local/lib
CFLAGS = -Wall $(DEBUG) $(GSLINC) -O3 -ffast-math -funroll-loops -I$(HOME)/.local/include
LFLAGS = -Wall $(DEBUG) $(GSL) -L$(HOME)/.local/lib/ -l:libgsl.so.28 -lgsl -lgslcblas -lm


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
	rm -fr *.o assemble
