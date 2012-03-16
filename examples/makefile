CPP = mpicxx
CC = mpicc
F90 = mpif90

INCLUDE = -I/u/nferraro/include -I$(HDF5_HOME)/include
LIBS = -L/u/nferraro/lib/_sunfire -lm3dc1 -L$(HDF5_HOME)/lib -lhdf5 -lz

EXAMPLES = cpp_example

all : $(EXAMPLES)

cpp_example : cpp_example.o
	$(CPP) $< $(LIBS) -o $@

clean :
	rm -fr *.o *.a *~ $(EXAMPLES)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@
