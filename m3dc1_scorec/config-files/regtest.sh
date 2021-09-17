# romulus
# module load gcc/4.8.5-v5m6xwi mpich/3.2.1-geowaxe cmake gsl/2.5-px4dg7h
# module unload zlib/1.2.11-vhzh5cf
#export ARCH=_romulus
#export BINDIR=/lore/seol/develop/M3DC1/unstructured
#export TESTDIR=$BINDIR/regtest
#export MPIRUN=mpirun

#stellar
# module load intel/2021.1 intel-mpi/intel/2021.1.1
# module load fftw/intel-19.1/intel-mpi/3.3.9
# module load hdf5/intel-2021.1/intel-mpi/1.10.6 gsl/2.6
#export ARCH=_stellar
#export BINDIR=/projects/M3DC1/scorec/src/M3DC1-trunk/unstructured
#export TESTDIR=/projects/M3DC1/scorec/src/M3DC1-trunk/unstructured/regtest
#export MPIRUN=srun
# salloc -n 64 -t 01:00:00

# centos7
# module load intel/2019.u3 openmpi/4.0.3 hdf5-parallel/1.10.5 fftw cmake git 
# module load simmodsuite/16.0-210626
# setenv ARCH _centos7
# setenv BINDIR /u/sseol/develop/M3DC1-trunk/unstructured
# setenv TESTDIR $BINDIR/regtest
# setenv MPIRUN srun
# salloc -n 64 -t 00:30:00 -p general --mem-per-cpu=2000M
#mpirun -n 16 /p/tsc/m3dc1/lib/SCORECLib/rhel7/intel2019u3-openmpi4.0.3/petsc3.13.5/bin/split_smb 

# sunfire
# module load intel/2019.u3 openmpi/4.0.3 hdf5-parallel/1.10.5 fftw cmake git 
# setenv ARCH _sunfire
# setenv BINDIR /u/sseol/develop/M3DC1-trunk/unstructured
# setenv BINDIR /u/sseol/develop/M3DC1-aug/unstructured
# setenv TESTDIR $BINDIR/regtest
# setenv MPIRUN srun
# salloc -n 64 -t 00:30:00 -p m3dc1 --mem-per-cpu=1000M
#mpirun -n 16 /p/tsc/m3dc1/lib/SCORECLib/rhel7/intel2019u3-openmpi4.0.3/petsc3.13.5/bin/split_smb 


cd  $TESTDIR/adapt/base
cp ../mesh/part* .
cp C1ke C1ke-org
$MPIRUN -n 16 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/KPRAD_2D/base
cp ../mesh/part* .
cp C1ke C1ke-org
$MPIRUN -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/KPRAD_restart/base
cp ../mesh/part* .
cp C1input.1 C1input
cp C1ke C1ke-org
$MPIRUN -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
cp C1input.2 C1input
$MPIRUN -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
rm -rf C1ke-org part*smb C1input

cd $TESTDIR/pellet/base
cp ../mesh/part* .
cp C1ke C1ke-org
$MPIRUN -n 64 $BINDIR/$ARCH-3d-opt-60/m3dc1_3d -ipetsc -options_file options_bjacobi.type_superludist
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/RMP/base
cp ../mesh/part* .
cp C1ke C1ke-org
$MPIRUN -n 16 $BINDIR/$ARCH-complex-opt-25/m3dc1_2d_complex -pc_factor_mat_solver_type mumps
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/RMP_nonlin/base
cp ../mesh/part* .
cp C1ke C1ke-org
$MPIRUN -n 64 $BINDIR/$ARCH-3d-opt-60/m3dc1_3d -ipetsc -options_file options_bjacobi.type_superludist
diff C1ke C1ke-org
rm -rf C1ke-org part*smb
