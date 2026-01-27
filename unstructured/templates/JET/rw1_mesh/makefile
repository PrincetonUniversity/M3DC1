all : input
	m3dc1_mfmgen input

fetch :
	wget -O jet_rw1.txt 'https://drive.google.com/uc?export=download&id=1YHaI5do8lMnFZVwyC8-RdJ5uSBae4NHQ'
	wget -O jet_rw1-15K0.smb 'https://drive.google.com/uc?export=download&id=1u7XyKh0OuD_6xFY9Mh4DzKh92N5IBjCM'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
