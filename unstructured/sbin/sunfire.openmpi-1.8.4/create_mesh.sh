#!/bin/bash

if [ ! -f input ]; then
    echo "Error: input file not found."
    exit
fi

#if command -v m3dc1_meshgen 2>/dev/null; then
#    echo "m3dc1_meshgen is present."

    m3dc1_meshgen input
#else
#    echo "m3dc1_meshgen is not present.  Trying create_smb."
#    POINTNAME=`awk '{ if($1 == "modelName")  print $2 }' input | tr '[:upper:]' '[:lower:]'`
#    THICKNESS=`awk '{ if($1 == "thickness")  print $2 }' input`
#    NAME=$POINTNAME$THICKNESS
#
#    echo $NAME
#
#    create_smb
#    convert_sim_sms $NAME.smd $NAME.sms $NAME.smb
#fi

