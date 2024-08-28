#!/bin/bash

MESH_FILE=mesh.smb
MODEL_FILE=model.dmg

if [ -z $1 ]; then
    echo "Usage: part_mesh.sh <model.smd> <mesh.sms> <parts>"
    echo "   or: part_mesh.sh <mesh.smb> <parts>"
else
    if [ ${1: -4} == ".smd" ]; then
	echo "Splitting .sms file"
	echo "*******************"
	echo "mesh convesion utility not available on Iris"
#	aprun -n 1 convert_sim_sms $1 $2 $MESH_FILE
#	aprun -n 1 make_model $MESH_FILE $MODEL_FILE
#	aprun -n $3 split_smb $MODEL_FILE $MESH_FILE part.smb $3
    elif [ ${1: -4} == ".smb" ]; then
	echo "Splitting .smb file"
	echo "*******************"
	mpiexec -n 1 make_model $1 $MODEL_FILE
	mpiexec -n $2 split_smb $MODEL_FILE $1 part.smb $2
    fi
fi

