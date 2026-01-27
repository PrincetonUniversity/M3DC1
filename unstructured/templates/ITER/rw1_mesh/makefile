all : input
	m3dc1_mfmgen input

fetch :
	wget -O iter_rw1.txt 'https://drive.google.com/uc?export=download&id=1qW2mUlljUSxhgWW91Nz4XGVjhrQBWfbU'
	wget -O iter_rw1-23K0.smb 'https://drive.google.com/uc?export=download&id=1CgHDjWAtFwzt5RUSO6fSvErgFaLBiWaj'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
