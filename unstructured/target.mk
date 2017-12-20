.SUFFIXES:

ifndef ARCH
  ARCH = $(shell echo $(HOST) | awk '{ sub(/[0-9]+/,""); print }')
# ARCH = $(shell uname -s)-$(shell uname -p)
endif

BIN_POSTFIX := $(ARCH)$(BIN_POSTFIX)

ifndef MAKECMDGOALS
 OBJDIR := _$(BIN_POSTFIX)
else
 OBJDIR := _$(ARCH)
endif

export BUILD_INFO = $(svn info)

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

.PHONY: cleanall
cleanall : 
	rm -fr _$(ARCH)*

.PHONY: clean
clean : 
	rm -fr _$(BIN_POSTFIX)*
