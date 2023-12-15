#!/bin/bash
#SBATCH --partition=all
#SBATCH --mem=128G
#SBATCH -n 48
#SBATCH -J M3DC1_regtest_KPRAD_2D
#SBATCH -t 00:30:00
#SBATCH -o C1stdout

touch started

source /opt/pppl/etc/profile.d/01-modules.sh

#create_fixed_mesh.sh AnalyticModel 0.06
#echo $? > mesh_created

PARTS=48
$M3DC1_MPIRUN -n $PARTS split_smb analytic-2K.smb part.smb $PARTS
echo $? > mesh_partitioned

$M3DC1_MPIRUN -n 48 m3dc1_2d #-pc_factor_mat_solver_type mumps -mat_mumps_icntl_14 100

touch finished
