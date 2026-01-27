all : input
	m3dc1_mfmgen input

fetch :
	wget -O east_fw1.txt 'https://drive.google.com/uc?export=download&id=17RJ7nanUP1F2_Y9XW9Vsdk-Mkfymer4x'
	wget -O east_fw1-12K0.smb 'https://drive.google.com/uc?export=download&id=1nynxKlrKJdDL9MaGW0VbgiuyUASTk8YF'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
