#!/bin/bash

if [ ! -e seed0.smb ]; then
    cp $M3DC1_DIR/bin/seed0.smb .
fi

srun -n 1 create_smb $1 $2
