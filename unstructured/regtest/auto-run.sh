#!/bin/csh

if ( $#argv == 0 ) then
    echo "Usage: auto-run.sh config-[arch].sh"
    exit(1)
endif

echo "================================================="
date
echo " _______________________"
hostname

source $argv[1]

cd $M3DC1_CODE_DIR
rm unstructured/regtest/run.log

git pull
cd unstructured
make OPT=1 clean; make OPT=1 COM=1 clean; make OPT=1 3D=1 MAX_PTS=60 clean
make OPT=1 
make OPT=1 COM=1 
make OPT=1 3D=1 MAX_PTS=60
make bin

cd regtest
./run $BATCH_SUFFIX > run.log
cat run.log
