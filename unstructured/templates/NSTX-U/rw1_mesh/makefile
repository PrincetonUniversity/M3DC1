all : input
	m3dc1_mfmgen input

.PHONY : clean all fetch

fetch :
	wget -O nstx-u_rw1.txt 'https://drive.google.com/uc?export=download&id=1hfol0kef8A453MXjzaoTjzdB9Pmx1z1s'
	wget -O nstx-u_rw1-9K0.smb 'https://drive.google.com/uc?export=download&id=1BeeauehEsOtHztw0MNCYfCO3obwqScL6'

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
