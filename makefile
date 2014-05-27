dirs = m3dc1_lib fusion_io trace_lib
#dirs = m3dc1_lib fusion_io
alldirs = $(dirs) examples

.PHONY : install clean $(alldirs)

all : $(dirs)

install : $(dirs)

shared : $(dirs)

clean : $(alldirs)
	rm -f *~

$(alldirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
#	mkdir -p $@/_$(M3DC1_ARCH)
#	$(MAKE) -C $@/_$(M3DC1_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
