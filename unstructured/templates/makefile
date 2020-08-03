dirs = AUG DIII-D NSTX NSTX-U KSTAR ITER JET

.PHONY : $(dirs) clean all

all : $(dirs)

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
