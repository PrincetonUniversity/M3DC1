.SUFFIXES:

OBJDIR := _$(M3DC1_ARCH)

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/makefile \
	SRCDIR=$(CURDIR) $(MAKECMDGOALS) 

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: clean
clean : 
	rm -fr *~ $(OBJDIR)
