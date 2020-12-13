BIN = m3dc1_skeleton
BIN_OMP = m3dc1_skeleton_omp
BIN_ACC = m3dc1_skeleton_acc
BIN_CUDA1 = m3dc1_skeleton_cuda1
BIN_CUDA = m3dc1_skeleton_cuda
BIN_OLD = m3dc1_skeleton_old

F90 = mpif90
LD = mpif90

ifeq ($(USEMPI), 1)
     #use xl compiler and netlib-lapack
     F90FLAGS = -qrealsize=8 -w -O3
     LDFLAGS := $(LDFLAGS)
endif
ifeq ($(USEOPENMP), 1)
     #use xl compiler and netlib-lapack
     F90FLAGS = -qrealsize=8 -w -O3
     F90FLAGS := $(F90FLAGS) -qsmp=omp -qoffload
     LDFLAGS := $(LDFLAGS) -qsmp=omp -qoffload
endif

ifeq ($(USEOPENACC), 1)
     #use pgi compiler
     F90FLAGS = -r8 -w -O3 -fast -Mlarge_arrays
     FFLAGS = -r8 -w -O3 -fast -Mlarge_arrays
     LDFLAGS = 

     #F90FLAGS := $(F90FLAGS) -acc -Minfo=accel -ta=multicore  #cc70,cuda10.0,time
     #FFLAGS := $(FFLAGS) -Mfixed -acc -Minfo=accel -ta=multicore  #cc70,cuda10.0,time
     F90FLAGS := $(F90FLAGS) -acc -Minfo=accel -ta=tesla:managed  #cc70,cuda10.0,time
     FFLAGS := $(FFLAGS) -Mfixed -acc -Minfo=accel -ta=tesla:managed  #cc70,cuda10.0,time
     LDFLAGS := $(LDFLAGS) -acc #-Mcuda
endif
ifeq ($(USECUDA), 1)
     #use pgi compiler
     F90FLAGS = -r8 -w -O3 -fast
     LDFLAGS = 
     F90FLAGS := $(F90FLAGS) -DUSECUDA
     LDFLAGS :=  $(LDFLAGS) -Mcuda #-L${MPI_ROOT}/lib -lmpi_ibm_mpifh
endif

LIB = -lblas -llapack
LIB_OMP = -L${OLCF_NETLIB_LAPACK_ROOT}/lib64 -lblas -llapack
LIB_CUDA = -lcudaforblas -lblas -llapack
INCLUDE = 

DEFS = -DUSE3D -DMAX_PTS=60

OBJS = etimes.o math.o element.o mesh.o nintegrate.o field.o \
	m3dc1_data.o m3dc1_nint.o matdef.o main.o

OBJS_OMP = etimes.o math.o element.o mesh.o nintegrate.o field.o \
	m3dc1_data.o m3dc1_nint.o matdef.o main_omp.o

OBJS_ACC = etimes.o math.o element_acc.o mesh_acc.o field_acc.o m3dc1_data_acc.o nintegrate_acc.o \
	dgemm.o dgemv.o m3dc1_nint_acc.o matdef_acc.o main_acc.o

OBJS_CUDA1 = etimes.o math.o element.o mesh.o field.o m3dc1_data.o nintegrate_cuda1.o \
	m3dc1_nint.o matdef_cuda1.o main_cuda.o

OBJS_CUDA = etimes.o math.o element.o mesh.o field.o m3dc1_data_cuda.o nintegrate_cuda.o \
	m3dc1_nint_cuda.o matdef_cuda.o main_cuda.o

OBJS_OLD = etimes.o math.o element.o mesh.o nintegrate_old.o field.o \
	m3dc1_data.o m3dc1_nint.o matdef_old.o main.o

.PHONY : clean dist

ifeq ($(USEOPENMP), 1)
all:	$(BIN_OMP)
endif
ifeq ($(USEMPI), 1)
all:	$(BIN)
endif
ifeq ($(USEOPENACC), 1)
all:	$(BIN_ACC)
endif
ifeq ($(USECUDA), 1)
#all:	$(BIN_CUDA1)
all:	$(BIN_CUDA)
endif

$(BIN) : $(OBJS)
	$(LD) $(LDFLAGS) $(LIB) $(OBJS) -o $@

$(BIN_OMP) : $(OBJS_OMP)
	$(LD) $(LDFLAGS) $(LIB_OMP) $(OBJS_OMP) -o $@

$(BIN_OLD) : $(OBJS_OLD)
	$(LD) $(LDFLAGS) $(LIB) $(OBJS_OLD) -o $@

$(BIN_ACC) : $(OBJS_ACC)
	$(LD) $(LDFLAGS) $(LIB) $(OBJS_ACC) -o $@

$(BIN_CUDA1) : $(OBJS_CUDA1)
	$(LD) $(LDFLAGS) $(LIB) $(OBJS_CUDA1) -o $@

$(BIN_CUDA) : $(OBJS_CUDA)
	$(LD) $(LDFLAGS) $(LIB_CUDA) $(OBJS_CUDA) -o $@

dist :
	tar c *.F90 README makefile > m3dc1_skeleton.tar
	gzip m3dc1_skeleton.tar

clean : 
	rm -fr *~ *.o *.mod
#rm -fr *~ *.o *.mod $(BIN) $(BIN_OLD) $(BIN_ACC) $(BIN_CUDA1) $(BIN_CUDA)

%.o : %.F90
	$(F90) -c $(INCLUDE) $(F90FLAGS) $(DEFS) $< -o $@

%.o : %.CUF
	$(F90) -c $(INCLUDE) $(F90FLAGS) $(DEFS) $< -o $@

%.o : %.f
	$(F90) -c $(INCLUDE) $(FFLAGS) $(DEFS) $< -o $@

