dirs = fw1_adapt rw1_adapt rw1_response rw1_stability

.PHONY : $(dirs) clean all

all : $(dirs)
	cp rw1_adapt/nstxu0.02.txt rw1_response
	cp rw1_adapt/nstxu0.02.txt rw1_stability

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
