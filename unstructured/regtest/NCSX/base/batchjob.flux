#!/bin/bash
#SBATCH --partition=all
#SBATCH --mem=128G
#SBATCH -n 64
#SBATCH -J regtest_NCSX
#SBATCH -t 00:30:00
#SBATCH -o C1stdout

touch started

source /opt/pppl/etc/profile.d/01-modules.sh

# partition mesh
PARTS=8
$M3DC1_MPIRUN -n $PARTS split_smb circle-858.smb part.smb $PARTS
echo $? > mesh_partitioned

# Run M3D-C1
$M3DC1_MPIRUN -n 64 m3dc1_3d_st -ipetsc -options_file options_bjacobi1

touch finished



