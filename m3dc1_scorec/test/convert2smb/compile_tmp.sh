#!/bin/sh -ex
export PKG_CONFIG_PATH=/fasttmp/fzhang/scorecbuild_new/install/lib/pkgconfig/
mpicxx main.cc `pkg-config --cflags --libs libapf_pumi`  -o main
#CC main.cc -I/global/project/projectdirs/mp288/edison/SCOREC-install/intel/include -Wl,--start-group -lapf -lapf_pumi -lpumi_util -lpumi_geom -lpcu -lpumi_geom_meshmodel -lpumi_mesh -lmeshadapt -Wl,--end-group -lzoltan -lparmetis -lmetis -L/global/project/projectdirs/mp288/edison/SCOREC-install/intel/lib -L/global/homes/h/hzhang/petsc/arch-xc30-opt/lib

