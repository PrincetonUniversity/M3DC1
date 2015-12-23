dirs = m3dc1_lib fusion_io trace_lib
alldirs = $(dirs) examples

.PHONY : install clean $(alldirs)

all : $(dirs)

install : $(dirs)

fio.tar.gz :
	tar c m3dc1_lib/*.h m3dc1_lib/*.cpp m3dc1_lib/makefile fusion_io/*.cpp fusion_io/*.h fusion_io/*.f90 fusion_io/*.F90 fusion_io/makefile install/* makefile > fio.tar
	gzip fio.tar

shared : $(dirs)

clean : $(alldirs)
	rm -f *~

$(alldirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
#	mkdir -p $@/_$(M3DC1_ARCH)
#	$(MAKE) -C $@/_$(M3DC1_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
