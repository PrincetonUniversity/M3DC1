dirs := $(wildcard */.)

.PHONY : $(dirs) clean all fetch

all : $(dirs)

fetch : $(dirs)

$(dirs) : 
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~
