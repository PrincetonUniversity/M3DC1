#!/bin/bash
#SBATCH --partition=all
#SBATCH --mem=128G
#SBATCH -n 16
#SBATCH -J M3DC1_regtest_RMP
#SBATCH -t 0:10:00
#SBATCH -o C1stdout

touch started

source /opt/pppl/etc/profile.d/01-modules.sh

#create_mesh.sh
#echo $? > mesh_created

PARTS=$SLURM_NTASKS
$M3DC1_MPIRUN -n $PARTS split_smb diiid-0.02-2.5-4.0-4K.smb part.smb $PARTS
echo $? > mesh_partitioned

$M3DC1_MPIRUN -n $SLURM_NTASKS m3dc1_2d_complex -pc_factor_mat_solver_type mumps

touch finished
