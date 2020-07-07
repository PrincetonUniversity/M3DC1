#!/bin/csh

echo "================================================="
date
echo " _______________________"
hostname

setenv M3DC1_CODE_DIR "/p/tsc/C1/src"
cd $M3DC1_CODE_DIR
rm unstructured/regtest/check.log

module purge
module use /p/m3dc1/modules
module load m3dc1/devel-centos7
module list

cd unstructured/regtest
./check centos7 > check.log
cat check.log

/usr/sbin/sendmail nferraro@pppl.gov < check.log
