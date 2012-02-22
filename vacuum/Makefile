 #=======================================================================
#
# Executable:	VACUUM
# Machine-OS:	i386-linux, alpha-osf1
# Author:	J.E. Menard, PPPL. Modified: M. Chance
# Date Written:	8/2002
#
#=======================================================================

TARGET	= vacuum

PARCOM	= 	\
#netcdf.inc	\
#vacdims.inc	\
vacuum1.inc	\
vacuum2.inc	\
vacuum3.inc	\
vacuum4.inc     \
vacuum5.inc	\
vacuum6.inc	\
vacuum7.inc	\
vacuum8.inc


OBJ77   = 	\
vacinc11.o      \
vacinc12.o      \
vacinc13.o      \
vacuum_ma.o   	\
vacuum_cal.o   	\
vacuum_diag.o 	\
vacuum_extl.o 	\
vacuum_hydra.o 	\
vacuum_io.o   	\
vacuum_math.o 	\
vacuum_pl.o   	\
vacuum_resh.o	\
vacuum_sens.o   \
vacuum_sprk2.o 	\
vacuum_ut.o   	\
vacuum_vac.o 

OBJTV   =       \
tv80lib.o

OBJ90	= 
OBJ95	= 
OBJCC	= bzio.o
#OBJCC	=

#=======================================================================
SHELL 	= /bin/csh

#vpath %.f /u/chance/Afs/Vacuum/Delta
#PTH = /u/chance/Afs/Vacuum/Delta
#vpath %.f /s/taurus/Chance/Delta
#PTH = /s/taurus/chance/Delta
vpath %.f /u/chance/Vacuum_SVN_Work/Vacuum_eps/trunk
PTH = /u/chance/Vacuum_SVN_Work/Vacuum_eps/trunk

#CCOMP	= cc
#FCOMP	= f77
#FCOPTS 	= 
#INCDIR 	= ./
#LDFLAGS = 
include makeflags.$(MACHTYPE).$(OSTYPE)
include makerules

#=======================================================================
