SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)

fetch : $(SUBDIRS)

clean: $(SUBDIRS)
	rm -f *~

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean fetch $(SUBDIRS)
