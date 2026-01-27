all : input
	m3dc1_mfmgen input

fetch :
	wget -O nstx_rw1.txt 'https://drive.google.com/uc?export=download&id=14VDRHuUGgpMKzGINpJNBkEFXXJcfAeln'
	wget -O nstx_rw1-10K0.smb 'https://drive.google.com/uc?export=download&id=1yyhw_6kYbFfdQUTnFelbxt5pHeRIF8NY'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
