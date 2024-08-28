$CXX main.cc -o main -I../../include -I../../api -Wl,-rpath,$SCORECROOT/lib -L$SCORECROOT/lib -Wl,-rpath,../../ -L../../  -L$SCORECROOT/lib -Wl,--start-group -lm3dc1_scorec $SCORECLIB -Wl,--end-group -L/fasttmp/fzhang/zoltan/lib -lzoltan -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib -lpetsc $PETSC_EXTERNAL_LIB_BASIC 

