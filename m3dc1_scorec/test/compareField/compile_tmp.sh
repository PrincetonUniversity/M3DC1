#!/bin/bash -ex
export SCORECROOT=/lore/fzhang/scorec_install/
export SCORECLIB="-lapf -lgmi -lm3dc1_scorec -lma -lparma  -lph -lapf_zoltan -lmds  -lpcu -lspr"
$CXX  main.cc -o main -I../../include  -I../../api -Wl,-rpath,$SCORECROOT/lib -L$SCORECROOT/lib -Wl,-rpath,../../.libs/ -L../../.libs/ -lm3dc1_scorec -L$SCORECROOT/lib -Wl,--start-group  $SCORECLIB -Wl,--end-group -L/fasttmp/fzhang/zoltan/lib -lzoltan -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib -lpetsc $PETSC_EXTERNAL_LIB_BASIC 
#./main tilt.txt mesh0.smb mesh0.smb field12 field12 1 1
