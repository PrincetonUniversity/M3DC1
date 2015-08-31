dirs = fw1 rw1

.PHONY : $(dirs) clean

$(dirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~