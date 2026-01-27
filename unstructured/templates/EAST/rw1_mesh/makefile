all : input
	m3dc1_mfmgen input

fetch :
	wget -O east_rw1.txt 'https://drive.google.com/uc?export=download&id=1o1-F2CJ-v2UjKdifFlOosFGer2AFH-NH'
	wget -O east_rw1-6K0.smb 'https://drive.google.com/uc?export=download&id=16adKaiXiC3s3Xfhuo0uc_xhh43oBtwBm'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
