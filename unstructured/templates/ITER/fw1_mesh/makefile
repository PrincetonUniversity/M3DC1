all : input
	m3dc1_mfmgen input

fetch :
	wget -O iter_fw1.txt 'https://drive.google.com/uc?export=download&id=1W9N61rnUSkfEyGnZ0mpRh01T0nwmyTz-'
	wget -O iter_fw1-11K0.smb 'https://drive.google.com/uc?export=download&id=1CXNJXuH5QwM-yorN46LzSHS2NeZP-taS'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
