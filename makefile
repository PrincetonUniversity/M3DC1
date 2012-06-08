dirs = m3dc1_lib fusion_io trace_lib

.PHONY : install clean $(dirs)

all : $(dirs)

install : $(dirs)

clean : $(dirs)

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
#	mkdir -p $@/_$(M3DC1_ARCH)
#	$(MAKE) -C $@/_$(M3DC1_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
