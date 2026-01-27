all : input
	m3dc1_mfmgen input

fetch :
	wget -O aug_fw1.txt 'https://drive.google.com/uc?export=download&id=1fTug1pMpbDBVpOMfTL3xQqu4XCtsy273'
	wget -O aug_fw1-9K0.smb 'https://drive.google.com/uc?export=download&id=1BXAE9nZnjxoUQL2ebyyv9japijLhZJk3'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
