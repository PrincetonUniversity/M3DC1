dirs = fw1_adapt fw1_response fw1_stability \
 rw1_adapt rw1_response rw1_stability

.PHONY : $(dirs) clean all

all : $(dirs)
	cp fw1_adapt/aug.txt fw1_response
	cp fw1_adapt/aug.txt fw1_stability
	cp rw1_adapt/aug-0.02-1.7-2.2.txt rw1_response
	cp rw1_adapt/aug-0.02-1.7-2.2.txt rw1_stability

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
