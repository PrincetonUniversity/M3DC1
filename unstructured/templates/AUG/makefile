dirs = fw1_mesh rw1_mesh

.PHONY : $(dirs) clean all fetch

all : $(dirs)

fetch : $(dirs)

$(dirs) :
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
