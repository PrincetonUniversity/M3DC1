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
install : install_idl install_doc
	echo $(ARCH)
	mkdir -m 755 -p $(INSTALL_DIR)
	mkdir -m 755 -p $(INSTALL_DIR)/batch
	-cp sbin/$(M3DC1_ARCH)/batch_script.* $(INSTALL_DIR)/batch
	-chmod 644 $(INSTALL_DIR)/batch/batch_script.* 
	mkdir -m 755 -p $(INSTALL_DIR)/bin
	cp sbin/extract_profiles.sh $(INSTALL_DIR)/bin
	chmod 755 $(INSTALL_DIR)/bin/extract_profiles.sh
	-cp sbin/$(M3DC1_ARCH)/*.sh $(INSTALL_DIR)/bin
	-chmod 755 $(INSTALL_DIR)/bin/*.sh
	-cp _$(ARCH)-opt-25/m3dc1_2d $(INSTALL_DIR)/bin
	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_2d
	-cp _$(ARCH)-complex-opt-25/m3dc1_2d_complex $(INSTALL_DIR)/bin
	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_2d_complex
	-cp _$(ARCH)-3d-opt-60/m3dc1_3d $(INSTALL_DIR)/bin
	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_3d

.PHONY: install_idl
install_idl : 
	mkdir -m 755 -p $(INSTALL_DIR)/idl
	cp idl/*.pro $(INSTALL_DIR)/idl
	chmod 644 $(INSTALL_DIR)/idl/*.pro

.PHONY: install_doc
install_doc :
	mkdir -m 755 -p $(INSTALL_DIR)/doc
	cp doc/* $(INSTALL_DIR)/doc
	-chmod 644 $(INSTALL_DIR)/doc/*

.PHONY: install_templates
install_templates : templates
	cp -r templates $(INSTALL_DIR)/
	find $(INSTALL_DIR) -type d -exec chmod 755 {} \;
	echo $(INSTALL_DIR)/templates/*/*_adapt | xargs -n 1 cp $(INSTALL_DIR)/batch/batch_script.adapt
	echo $(INSTALL_DIR)/templates/*/*_response $(INSTALL_DIR)/templates/*/*_stability | xargs -n 1 cp $(INSTALL_DIR)/batch/batch_script.2d_complex
	find $(INSTALL_DIR)/templates -type f -exec chmod 644 {} \;
