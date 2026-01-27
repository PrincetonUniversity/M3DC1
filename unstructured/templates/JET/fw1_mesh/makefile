all : input
	m3dc1_mfmgen input

fetch :
	wget -O jet_fw1.txt 'https://drive.google.com/uc?export=download&id=1cXIpZIHGh0IPoLHZcFSwz_TpWT6L-2BG'
	wget -O jet_fw1-11K0.smb 'https://drive.google.com/uc?export=download&id=136hXF6gGzHZtT1Ph0les8O9JRiD9s2ez'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
