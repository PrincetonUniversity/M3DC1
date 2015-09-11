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

VERSION = $(shell cat release_version)
INSTALL_DIR = $(M3DC1_INSTALL_DIR)/m3dc1-$(M3DC1_ARCH)-$(VERSION)

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
	cd templates ; make clean

.PHONY: clean
clean : 
	rm -fr _$(BIN_POSTFIX)*

.PHONY: templates
templates :
	cd templates; make

.PHONY: install
install : templates
	echo $(ARCH)
	mkdir -p $(INSTALL_DIR)
	mkdir -p $(INSTALL_DIR)/idl
	cp idl/*.pro $(INSTALL_DIR)/idl
	mkdir -p $(INSTALL_DIR)/bin
	cp sbin/extract_profiles.sh $(INSTALL_DIR)/bin
	cp sbin/part_mesh.$(M3DC1_ARCH).sh $(INSTALL_DIR)/bin/part_mesh.sh
	cp sbin/create_fixed_mesh.$(M3DC1_ARCH).sh $(INSTALL_DIR)/bin/create_fixed_mesh.sh
	cp sbin/create_mesh.$(M3DC1_ARCH).sh $(INSTALL_DIR)/bin/create_mesh.sh
	cp _$(ARCH)-opt-25/m3dc1_2d $(INSTALL_DIR)/bin
	cp _$(ARCH)-complex-opt-25/m3dc1_2d_complex $(INSTALL_DIR)/bin
	cp _$(ARCH)-3d-opt-60/m3dc1_3d $(INSTALL_DIR)/bin
	cp -r templates $(INSTALL_DIR)/
