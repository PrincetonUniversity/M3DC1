dirs = DIII-D NSTX-U

.PHONY : $(dirs) clean all

all : $(dirs)

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
