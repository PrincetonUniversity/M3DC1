#!/bin/bash

ERROR=1
if [ -e "$1" ]; then   
if [[ "$1" == m3dc1_profiles_*.txt ]]; then
    echo "Reading Shafer profiles."

    tail -n +2 $1 | awk '{print $1 " " $2}' > profile_te
    tail -n +2 $1 | awk '{print $1 " " $3*0.1}' > profile_ne
    tail -n +2 $1 | awk '{print $1 " " $4}' > profile_omega

    ERROR=0

elif [[ "$1" == *.tgz  ]]; then

    echo "Reading tarball from Osborne's tools."

    # create a directory to unpack profiles
    mkdir -p profiles

    # unpack the tarball 
    echo "Unpacking tarball"
    tar xfz $1 -C profiles/

    # write profiles into format readable by M3D-C1
    echo "Extracting profiles"
    tail -n +2 profiles/netanh*psi_*.dat > profile_ne
    tail -n +2 profiles/tetanh*psi_*.dat > profile_te
    tail -n +2 profiles/omgebspl*psi_*.dat > profile_omega

    ERROR=0

elif [[ "$1" == p*.* ]]; then
    echo "Reading p-eqdsk file"

    sed -n '/psinorm ne/,/psinorm/{/psinorm/!p}' $1 > profile_ne
    sed -n '/psinorm te/,/psinorm/{/psinorm/!p}' $1 > profile_te
    sed -n '/psinorm omgeb/,/psinorm/{/psinorm/!p}' $1 > profile_omega

    ERROR=0
fi
fi

if [ $ERROR != 0 ]; then
    echo "Usage: ./extract_profiles <profile_file>"
    echo "  where <profile_file> is one of: "
    echo "  * a m3dc1_profiles_*.txt from Shafer"
    echo "  * a .tgz file from Osborne's phython tools"
    echo "  * a p-eqdsk file"
    exit $ERROR
else
    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"
fi