#!/bin/bash

SCOREC_DIR=/global/project/projectdirs/mp288/hopper/scorec/utilities/create_smb
SEED_MESH=seed0.smb
CREATE_SMB=create_smb 

if [ ! -e $SEED_MESH ]; then
    cp $SCOREC_DIR/$SEED_MESH .
fi

aprun -n 1 $SCOREC_DIR/$CREATE_SMB $1 $2
