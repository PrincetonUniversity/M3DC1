dirs = nstx1 nstx2 nstx3 nstx4 kink1

ARCH = $(shell echo $(HOST) | awk '{ sub(/[0-9]+/,""); print }')
M3DC1 = mpiexec -n $(N) ../../_$(ARCH)
blue='\e[0;34m'
NC='\e[0m'

export M3DC1

.PHONY : clean $(dirs)

all : $(dirs)

clean : $(dirs)

reset : $(dirs)

$(dirs) : 
	@echo -e ${blue}"$@..."${NC}
	@cd $@; $(MAKE) $(MAKECMDGOALS)