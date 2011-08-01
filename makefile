dirs = m3dc1_lib trace_lib convert

.PHONY : install clean $(dirs)

all : $(dirs)

install : $(dirs)

clean : 
	cd m3dc1_lib ; make clean
	cd trace_lib ; make clean
	cd convert ; make clean

$(dirs) : 
	mkdir -p $@/_$(M3DC1_ARCH)
	$(MAKE) -C $@/_$(M3DC1_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
