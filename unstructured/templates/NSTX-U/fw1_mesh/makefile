all : input
	m3dc1_mfmgen input

fetch :
	wget -O nstx-u_fw1.txt 'https://drive.google.com/uc?export=download&id=1HRyK6fo-6PriF9v1ZhGUyCJCAOMiBbUL'
	wget -O nstx-u_fw1-14K0.smb 'https://drive.google.com/uc?export=download&id=1926_1_bIkaVniUtWnCDvzXF8kPK3NDwN'

.PHONY : clean all fetch
clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
