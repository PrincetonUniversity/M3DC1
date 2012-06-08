CC = mpicxx CC='icc'
CFLAGS = -O3 -Wall

F90 = mpif90
F90FLAGS = -r8 -O3 -Minform=warn

LD = $(F90)
LDFLAGS = 

ARCH = ar cr

HDF5_LIBS = -L/u/ferraro/hdf5-1.8.6-linux-i386/lib -lhdf5 -lz \
	-Wl,-rpath,/u/ferraro/hdf5-1.8.6-linux-i386/lib
HDF5_INCLUDE = -I/u/ferraro/hdf5-1.8.6-linux-i386/include

INSTALL_DIR = /u/ferraro/bin

LIBS = -lstdc++ -lpmpich++
