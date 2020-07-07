#!/bin/csh

echo "================================================="
date
echo " _______________________"
hostname

setenv M3DC1_CODE_DIR "/p/tsc/C1/src"
cd $M3DC1_CODE_DIR
rm unstructured/regtest/run.log

module purge
module use /p/m3dc1/modules
module load m3dc1/devel-centos7
module load git
module list
printenv PATH

setenv GIT_SSH_COMMAND 'ssh -i /u/nferraro/.ssh/id_rsa_0 -l nferraro' 
git pull
cd unstructured
make OPT=1 clean; make OPT=1 COM=1 clean; make OPT=1 3D=1 MAX_PTS=60 clean
make OPT=1 
make OPT=1 COM=1 
make OPT=1 3D=1 MAX_PTS=60
make bin

cd regtest
./run centos7 > run.log
cat run.log
