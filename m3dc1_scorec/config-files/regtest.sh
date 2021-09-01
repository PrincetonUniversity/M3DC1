salloc -n 64 -t 01:00:00

#stellar
# module load intel/2021.1 intel-mpi/intel/2021.1.1
# module load fftw/intel-19.1/intel-mpi/3.3.9
# module load hdf5/intel-2021.1/intel-mpi/1.10.6 gsl/2.6
#export ARCH=_stellar-seol
#export BINDIR=/projects/M3DC1/scorec/src/M3DC1-trunk/unstructured
#export TESTDIR=/projects/M3DC1/scorec/src/M3DC1-trunk/unstructured/regtest

cd  $TESTDIR/adapt/base
cp ../mesh/part* .
cp C1ke C1ke-org
srun -n 16 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/KPRAD_2D/base
cp ../mesh/part* .
cp C1ke C1ke-org
srun -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
srun -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
rm -rf C1ke-org part*smb

cd $TESTDIR/KPRAD_restart/base
cp ../mesh/part* .
cp C1input.1 C1input
cp C1ke C1ke-org
srun -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
cp C1input.2 C1input
srun -n 48 $BINDIR/$ARCH-opt-25/m3dc1_2d
diff C1ke C1ke-org
rm -rf C1ke-org part*smb C1input

cd $TESTDIR/pellet/base
cp ../mesh/part* .
cp C1ke C1ke-org
srun -n 64 $BINDIR/$ARCH-3d-opt-60/m3dc1_3d -ipetsc -options_file options_bjacobi.type_superludist
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/RMP/base
cp ../mesh/part* .
cp C1ke C1ke-org
srun -n 16 $BINDIR/$ARCH-complex-opt-25/m3dc1_2d_complex -pc_factor_mat_solver_type mumps
diff C1ke C1ke-org
rm -rf C1ke-org part*smb

cd $TESTDIR/RMP_nonlin/base
cp ../mesh/part* .
cp C1ke C1ke-org
srun -n 64 $BINDIR/$ARCH-3d-opt-60/m3dc1_3d -ipetsc -options_file options_bjacobi.type_superludist
diff C1ke C1ke-org
rm -rf C1ke-org part*smb
