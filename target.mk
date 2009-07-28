.SUFFIXES:

ifndef ARCH
 ARCH = $(shell uname -s)-$(shell uname -p)
endif

BIN_POSTFIX := $(ARCH)$(BIN_POSTFIX)

ifndef MAKECMDGOALS
 OBJDIR := _$(BIN_POSTFIX)
else
 OBJDIR := _$(ARCH)
endif

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/makefile \
	SRCDIR=$(CURDIR) ARCH=$(ARCH) BIN_POSTFIX=$(BIN_POSTFIX) \
	$(MAKECMDGOALS) 

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: clean
clean : 
	rm -fr _$(ARCH)* ../common/*.o
