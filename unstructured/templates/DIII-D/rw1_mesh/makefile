all : input
	m3dc1_mfmgen input

fetch :	
	wget -O diii-d_rw1.txt 'https://drive.google.com/uc?export=download&id=1D8w4xWb5GBPfGJE5QhIE7q9xkVteePqx'
	wget -O diii-d_rw1-9K0.smb 'https://drive.google.com/uc?export=download&id=1nouORw813fhOVij80apVHluErpYMVUXy'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
