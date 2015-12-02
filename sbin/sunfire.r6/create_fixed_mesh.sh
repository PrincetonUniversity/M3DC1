#!/bin/bash

SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/utilities/create_smb
SEED_MESH=seed0.smb
CREATE_SMB=create_smb 

if [ ! -e $SEED_MESH ]; then
    cp $SCOREC_DIR/$SEED_MESH .
fi

$SCOREC_DIR/$CREATE_SMB $1 $2
