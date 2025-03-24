#!/bin/bash

CONVERT=/p/tsc/m3dc1/lib/SCORECLib/rhel6/utilities/create_mesh/convert_sim_sms
SPLIT_SMB=/global/project/projectdirs/mp288/hopper/scorec/utilities/split_smb/split_smb
MAKE_MODEL=/global/project/projectdirs/mp288/hopper/scorec/utilities/split_smb/make_model
MESH_FILE=mesh.smb
MODEL_FILE=model.dmg

if [ -z $1 ]; then
    echo "Usage: part_mesh.sh <model.smd> <mesh.sms> <parts>"
    echo "   or: part_mesh.sh <mesh.smb> <parts>"
else
    if [ ${1: -4} == ".smd" ]; then
	echo "Splitting .sms file"
	echo "*******************"
	aprun -n 1 $CONVERT $1 $2 $MESH_FILE
	aprun -n 1 $MAKE_MODEL $MESH_FILE $MODEL_FILE
	aprun -n $3 $SPLIT_SMB $MODEL_FILE $MESH_FILE part.smb $3
    elif [ ${1: -4} == ".smb" ]; then
	echo "Splitting .smb file"
	echo "*******************"
	aprun -n 1 $MAKE_MODEL $1 $MODEL_FILE
	aprun -n $2 $SPLIT_SMB $MODEL_FILE $1 part.smb $2
    fi
fi

