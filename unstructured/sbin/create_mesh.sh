#!/bin/bash

meshgen=$(which m3dc1_meshgen)

if [ -x "$meshgen" ]; then
    echo "m3dc1_meshgen found."
    m3dc1_meshgen input
else
    echo "m3dc1_meshgen not found.  Trying create_smb"
    POINTNAME=`awk '{ if($1 == "modelName")  print $2 }' input | tr '[:upper:]' '[:lower:]'`
    THICKNESS=`awk '{ if($1 == "thickness")  print $2 }' input`
    NAME=$POINTNAME$THICKNESS
    echo $NAME

    create_smb
    convert_sim_sms $NAME.smd $NAME.sms $NAME.smb
fi
