.SUFFIXES:

ifndef ARCH
  ifdef M3DC1_ARCH
    ARCH = $(M3DC1_ARCH)
  else
    ARCH = $(shell echo $(HOST) | awk '{ sub(/[0-9]+/,""); print }')
  endif
endif

BIN_POSTFIX := $(ARCH)$(BIN_POSTFIX)

ifndef MAKECMDGOALS
 OBJDIR := _$(BIN_POSTFIX)
else
 OBJDIR := _$(ARCH)
endif

VERSION = $(shell cat release_version)
INSTALL_DIR = $(M3DC1_INSTALL_DIR)/m3dc1-$(ARCH)-$(VERSION)

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/makefile \
	SRCDIR=$(CURDIR) ARCH=$(ARCH) BIN_POSTFIX=$(BIN_POSTFIX) \
	$(MAKECMDGOALS) 

include $(ARCH).mk

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: all
all :
	make OPT=1
	make OPT=1 COM=1
	make OPT=1 COM=1 PAR=1
	make OPT=1 3D=1 MAX_PTS=60
	make OPT=1 3D=1 MAX_PTS=60 PAR=1
	make OPT=1 3D=1 MAX_PTS=125 ST=1
	make OPT=1 3D=1 MAX_PTS=125 ST=1 PAR=1
	make a2cc
	make bin

.PHONY: cleanall
cleanall : 
	rm -fr _$(ARCH)*
	cd templates ; make clean

.PHONY: clean
clean : 
	rm -fr _$(ARCH)
	rm -fr _$(ARCH)-*

.PHONY: templates
templates :
	cd templates; make

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

.PHONY: install_tutorials
install_tutorials :
	mkdir -m 755 -p $(INSTALL_DIR)/tutorials
	cp -r tutorials/* $(INSTALL_DIR)/tutorials
	find $(INSTALL_DIR)/tutorials -type d -exec chmod 755 {} \;
	find $(INSTALL_DIR)/tutorials -type f -exec chmod 644 {} \;

.PHONY: install_templates
install_templates : templates
	cp -r templates $(INSTALL_DIR)/
	find $(INSTALL_DIR) -type d -exec chmod 755 {} \;
	echo $(INSTALL_DIR)/templates/*/*_adapt | xargs -n 1 cp $(INSTALL_DIR)/batch/batch_script.adapt
	echo $(INSTALL_DIR)/templates/*/*_response $(INSTALL_DIR)/templates/*/*_stability | xargs -n 1 cp $(INSTALL_DIR)/batch/batch_script.2d_complex
	find $(INSTALL_DIR)/templates -type f -exec chmod 644 {} \;

.PHONY: install_device_data
install_device_data : # device_data
	mkdir -m 755 -p $(INSTALL_DIR)/device_data
	cp -r device_data/* $(INSTALL_DIR)/device_data
	find $(INSTALL_DIR)/device_data -type d -exec chmod 755 {} \;
	find $(INSTALL_DIR)/device_data -type f -exec chmod 644 {} \;

#.PHONY: install_scorec
#install_scorec : 
#	mkdir -m 755 -p $(INSTALL_DIR)/bin
#	echo $(SCOREC_UTIL_DIR)
#	-cp $(SCOREC_UTIL_DIR)/create_smb/create_smb $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/create_smb
#	-cp $(SCOREC_UTIL_DIR)/seed0.smb $(INSTALL_DIR)/bin
#	-chmod 644 $(INSTALL_DIR)/bin/seed0.smb
#	cp $(SCOREC_UTIL_DIR)/split_smb/split_smb $(INSTALL_DIR)/bin
#	chmod 755 $(INSTALL_DIR)/bin/split_smb
#	cp $(SCOREC_UTIL_DIR)/split_smb/make_model $(INSTALL_DIR)/bin
#	chmod 755 $(INSTALL_DIR)/bin/make_model
#	cp $(SCOREC_UTIL_DIR)/split_smb/make_model $(INSTALL_DIR)/bin
#	chmod 755 $(INSTALL_DIR)/bin/make_model
#	-cp $(SCOREC_UTIL_DIR)/m3dc1_meshgen/m3dc1_meshgen $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_meshgen
#	-cp $(SCOREC_UTIL_DIR)/m3dc1_meshgen/convert_sim_sms $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/convert_sim_sms

.PHONY: install
install : install_idl install_doc install_device_data install_tutorials #install_scorec
	echo $(ARCH)
	mkdir -m 755 -p $(INSTALL_DIR)
	mkdir -m 755 -p $(INSTALL_DIR)/batch
	-cp sbin/$(ARCH)/batch_script.* $(INSTALL_DIR)/batch
	-chmod 644 $(INSTALL_DIR)/batch/batch_script.* 
	mkdir -m 755 -p $(INSTALL_DIR)/bin
	-cp _$(ARCH)/bin/* $(INSTALL_DIR)/bin
	chmod 755 $(INSTALL_DIR)/bin/*
#	cp sbin/extract_profiles.sh $(INSTALL_DIR)/bin
#	chmod 755 $(INSTALL_DIR)/bin/extract_profiles.sh
#	cp sbin/m3dc1_units.sh $(INSTALL_DIR)/bin
#	chmod 755 $(INSTALL_DIR)/bin/m3dc1_units.sh
#	-cp sbin/$(ARCH)/*.sh $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/*.sh
#	-cp _$(ARCH)/a2cc $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/a2cc
#	-cp _$(ARCH)-opt-25/m3dc1_2d $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_2d
#	-cp _$(ARCH)-complex-opt-25/m3dc1_2d_complex $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_2d_complex
#	-cp _$(ARCH)-3d-opt-60/m3dc1_3d $(INSTALL_DIR)/bin
#	-chmod 755 $(INSTALL_DIR)/bin/m3dc1_3d
