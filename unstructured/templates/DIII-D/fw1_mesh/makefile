all : input
	m3dc1_mfmgen input

fetch :
	wget -O diii-d_fw1.txt 'https://drive.google.com/uc?export=download&id=136kEaygMe-ACKGyAYrVShgvMZ7IAAGMY'
	wget -O diii-d_fw1-15K0.smb 'https://drive.google.com/uc?export=download&id=1q0LQ_tM4fEYJZwFuuwfzJQCXH1Ja5G8k'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
