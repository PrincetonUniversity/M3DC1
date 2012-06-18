dirs = m3dc1_lib fusion_io trace_lib
alldirs = $(dirs) examples

.PHONY : install clean $(alldirs)

all : $(dirs)

install : $(dirs)

clean : $(alldirs)

$(alldirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
#	mkdir -p $@/_$(M3DC1_ARCH)
#	$(MAKE) -C $@/_$(M3DC1_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
