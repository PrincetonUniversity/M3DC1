#!/bin/bash

POINTNAME=`awk '{ if($1 == "modelName")  print $2 }' input | tr '[:upper:]' '[:lower:]'`
THICKNESS=`awk '{ if($1 == "thickness")  print $2 }' input`
NAME=$POINTNAME$THICKNESS

echo $NAME

/home/jinchen/lib/scorec/utilities/create_mesh/create_smd
/home/jinchen/lib/scorec/utilities/create_mesh/convert_sim_sms $NAME.smd $NAME.sms $NAME.smb
