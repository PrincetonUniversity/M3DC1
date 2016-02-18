dirs = AUG DIII-D NSTX NSTX-U KSTAR

.PHONY : $(dirs) clean all

all : $(dirs)

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
