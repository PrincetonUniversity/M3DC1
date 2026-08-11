#!/bin/bash
#SBATCH --partition=all
#SBATCH --mem=128G
#SBATCH -N 1
#SBATCH -J m3dc1_regtest_adapt
#SBATCH -t 0:10:00
#SBATCH -o C1stdout

touch started

source /opt/pppl/etc/profile.d/01-modules.sh

# partition mesh
PARTS=16
$M3DC1_MPIRUN -n $PARTS split_smb circle-0.02-0.0-0.0-2K.smb part.smb $PARTS
echo $? > mesh_partitioned

# Run M3D-C1
$M3DC1_MPIRUN -n 16 m3dc1_2d
$M3DC1_MPIRUN -n 16 collapse circle-0.02-0.0-0.0.dmg ts0-adapted.smb mesh.smb $SLURM_NTASKS

if ! [ -f mesh0.smb ]; then
    echo 0 > C1ke
fi

touch finished
