#!/bin/bash


if [ -e "$1" ]; then   
if [[ "$1" == m3dc1_profiles_*.txt ]]; then
    echo "Reading Shafer profiles."

    tail -n +2 $1 | awk '{print $1 " " $2}' > profile_te
    tail -n +2 $1 | awk '{print $1 " " $3*0.1}' > profile_ne
    tail -n +2 $1 | awk '{print $1 " " $4}' > profile_omega

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"

    exit 0

elif [[ "$1" == *_ipec_profs.dat ]]; then
    echo "Reading IPEC profiles."

    tail -n +7 $1 | awk '{print $1 " " $5*1e-3}' > profile_te
    tail -n +7 $1 | awk '{print $1 " " $3*1e-20}' > profile_ne
    tail -n +7 $1 | awk '{print $1 " " $6}' > profile_omega

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"

    exit 0

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
    tail -n +2 profiles/omgebspl*psi_*.dat > profile_omega.ExB
    tail -n +2 profiles/ommvbspl*psi_*.dat > profile_omega.ion
    cp profile_omega.ExB profile_omega
    
    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"

    exit 0

elif [[ "$1" == *_ntvin.dat ]]; then
    echo "Reading NTVIN file"

    tail -n +11 $1 | awk '{print $2 " " $3*1e-3}' > profile_omega
    tail -n +11 $1 | awk '{print $2 " " $5*1e-20}' > profile_ne
    tail -n +11 $1 | awk '{print $2 " " $9*1e-3}' > profile_te

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega = 1"

    exit 0

elif [[ "$1" == k.* ]]; then
    echo "Reading k. file"

    tail -n +11 $1 | awk '{print $1 " " $8*1e-3}' > profile_omega
    tail -n +11 $1 | awk '{print $1 " " $4*1e-20}' > profile_ne
    tail -n +11 $1 | awk '{print $1 " " $6*1e-3}' > profile_te

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"

    exit 0

elif [[ "$1" == NSTX*_profiles.dat ]]; then
    echo "Reading NSTX profiles file"

    tail -n +7 $1 | awk '{print $1*$1 " " $6}' > profile_omega
    tail -n +7 $1 | awk '{print $1*$1 " " $2*1e-20}' > profile_ne
    tail -n +7 $1 | awk '{print $1*$1 " " $3*1e-3}' > profile_te

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega = 1"

    exit 0

elif [[ "$1" == pdbne*.dat ]]; then
    echo "Reading pdbne file"

    tail -n +2 $1 | awk '{print $1 " " $2*1e-6}' > profile_ne_rho_0
    echo "In C1input set"
    echo " iread_ne = 4"
    exit 0

elif [[ "$1" == pdbte*.dat ]]; then
    echo "Reading pdbte file"

    tail -n +2 $1 | awk '{print $1 " " $2}' > profile_te_rho_3
    echo "In C1input set"
    echo " iread_te = 4"
    exit 0

elif [[ "$1" == pdbomgeb*.dat ]]; then
    echo "Reading pdbomgeb file"

    tail -n +2 $1 | awk '{print $1 " " $2*1e3}' > profile_omega_rho_0
    echo "In C1input set"
    echo " iread_omega_ExB = 4"
    exit 0

elif [[ "$1" == p*.* ]]; then
    echo "Reading p-eqdsk file"

    sed -n '/psinorm ne/,/psinorm/{/psinorm/!p}' $1 > profile_ne
    sed -n '/psinorm te/,/psinorm/{/psinorm/!p}' $1 > profile_te
    sed -n '/psinorm omgeb/,/psinorm/{/psinorm/!p}' $1 > profile_omega.ExB
    sed -n '/psinorm omgvb/,/psinorm/{/psinorm/!p}' $1 > profile_omega.C
    sed -n '/psinorm ommvb/,/psinorm/{/psinorm/!p}' $1 > profile_omega.ion
    sed -n '/psinorm omevb/,/psinorm/{/psinorm/!p}' $1 > profile_omega.electron
    cp profile_omega.ExB profile_omega

    echo "In C1input set"
    echo " iread_ne = 1"
    echo " iread_te = 1"
    echo " iread_omega_ExB = 1"

    exit 0

elif [[ "$1" == neprof_*.asc ]]; then
    echo "Reading neprof.asc file"

    sed 's/D/E/g' $1 | awk '{print $1*$1 " " $2*1e-20}' > profile_ne
    echo "In C1input set"
    echo " iread_ne = 1"
    exit 0

elif [[ "$1" == Teprof_*.asc ]]; then
    echo "Reading Teprof.asc file"

    sed 's/D/E/g' $1 | awk '{print $1*$1 " " $2*1e-3}' > profile_te
    echo "In C1input set"
    echo " iread_te = 1"
    exit 0

elif [[ "$1" == vtprof_*.asc ]]; then
    echo "Reading vtprof.asc file"

    sed 's/D/E/g' $1 | awk '{print $1*$1 " " $2}' > profile_vphi
    echo "In C1input set"
    echo " iread_omega = 3"
    exit 0

elif [[ "$1" == *_profiles.dat ]]; then
    echo "Reading AUG profile file"
    tail -n +2 $1 | awk '{printf "%10g\t%10g\n", $1*$1, $8*1e-4}' > profile_ne
    tail -n +2 $1 | awk '{printf "%10g\t%10g\n", $1*$1, $10*1e-3}' > profile_te
    echo "In C1input set"
    echo " iread_te = 1"
    echo " iread_ne = 1"
    exit

elif [[ "$1" == *_vtor_*.dat ]]; then
    echo "Reading AUG rotation profile file"
    tail -n +2 $1 | awk '{printf "%10g\t%10g\n", $3, $4*1e-3/$1}' > profile_omega
    echo "In C1input set"
    echo " iread_omega = 1"
    exit

fi
fi

echo "Usage: ./extract_profiles <profile_file>"
echo "  where <profile_file> is one of: "
echo "  * a m3dc1_profiles_*.txt from Shafer"
echo "  * a .tgz file from Osborne's phython tools"
echo "  * a p-eqdsk file"
echo "  * a ntvin.dat file"
echo "  * a NSTX_profiles.dat file"
echo "  * a pdbne*.dat, pdbte*.dat, or pdbomgeb*.dat file"

exit 1
