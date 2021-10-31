#!/bin/bash

MESH_FILE=mesh.smb
MODEL_FILE=model.dmg
RUN=$(which "$M3DC1_MPIRUN")

if [ -z $1 ]; then
    echo "Usage: part_mesh.sh <model.smd> <mesh.sms> <parts> <old_parts>"
    echo "   or: part_mesh.sh <mesh.smb> <parts> <old_parts>"
else
    if [ -x "$RUN" ]; then
	echo "found it: "

	if [ ${1: -4} == ".smd" ]; then
	    echo "Splitting .sms file"
	    echo "*******************"
	    $M3DC1_MPIRUN -n 1 convert_sim_sms $1 $2 $MESH_FILE
	    $M3DC1_MPIRUN -n 1 make_model $MESH_FILE $MODEL_FILE
	    $M3DC1_MPIRUN -n $3 split_smb $MODEL_FILE $MESH_FILE part.smb $3
	elif [ ${1: -4} == ".smb" ]; then
	    echo "Splitting .smb file"
	    echo "*******************"
	    $M3DC1_MPIRUN -n 1 make_model $1 $MODEL_FILE
	    $M3DC1_MPIRUN -n $2 split_smb $MODEL_FILE $1 part.smb $2
	fi
    else
	echo "ERROR: M3DC1_MPIRUN environment variable not set."
	echo "       This should be set to the command for running MPI jobs"
	echo "       (e.g. mpiexec or srun)"
    fi
fi

