#!/bin/bash
#SBATCH -n 16
#SBATCH -J M3DC1_Tutorial
#SBATCH -t 2:00:00
#SBATCH --mem=180000

PARTS=$SLURM_NTASKS
$M3DC1_MPIRUN -n $PARTS split_smb mesh.smb part.smb $PARTS

$M3DC1_MPIRUN -n $SLURM_NTASKS m3dc1_2d_complex -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 100
