dirs = fw1_adapt fw1_response fw1_stability \
 rw1_adapt rw1_response rw1_stability

.PHONY : $(dirs) clean all

all : $(dirs)
	cp fw1_adapt/jet.txt fw1_response
	cp fw1_adapt/jet.txt fw1_stability
	cp rw1_adapt/jet-0.02-2.5-4.0.txt rw1_response
	cp rw1_adapt/jet-0.02-2.5-4.0.txt rw1_stability

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
