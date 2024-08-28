#!/bin/bash

POINTNAME=`awk '{ if($1 == "modelName")  print $2 }' input | tr '[:upper:]' '[:lower:]'`
THICKNESS=`awk '{ if($1 == "thickness")  print $2 }' input`
NAME=$POINTNAME$THICKNESS

echo $NAME

echo "SCOREC create_mesh utility not available on SATURN"
#/global/project/projectdirs/mp288/edison/scorec/utilities/create_mesh/create_smb
#/global/project/projectdirs/mp288/edison/scorec/utilities/create_mesh/convert_sim_sms $NAME.smd $NAME.sms $NAME.smb
