#!/bin/bash
#SBATCH --partition=all
#SBATCH --mem=128G
#SBATCH -n 64
#SBATCH -J M3DC1_regtest_RMP_nonlin
#SBATCH -t 0:30:00
#SBATCH -o C1stdout

touch started

source /opt/pppl/etc/profile.d/01-modules.sh

# partition mesh
PARTS=16
$M3DC1_MPIRUN -n $PARTS split_smb circle-0.10-0.0-0.0-1K.smb part.smb $PARTS
echo $? > mesh_partitioned

# Run M3D-C1
$M3DC1_MPIRUN -n 64 m3dc1_3d -ipetsc -options_file options_bjacobi.type_superludist

touch finished
