all : input
	m3dc1_mfmgen input

fetch :
	wget -O kstar_fw1.txt 'https://drive.google.com/uc?export=download&id=1G6pG_LkTWXg3wp59N_V0wL431H-OzS-O'
	wget -O kstar_fw1-10K0.smb 'https://drive.google.com/uc?export=download&id=1BL3UYqWhUaxotP5YCCUJSYfYhd2rEjpZ'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
