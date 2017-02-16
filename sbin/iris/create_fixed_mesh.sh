#!/bin/bash

SCOREC_DIR=/fusion/projects/codes/m3dc1/scorec/tools/create_smb
SEED_MESH=seed0.smb
CREATE_SMB=create_smb 

if [ ! -e $SEED_MESH ]; then
    cp $SCOREC_DIR/$SEED_MESH .
fi

mpiexec -n 1 $SCOREC_DIR/$CREATE_SMB $1 $2
