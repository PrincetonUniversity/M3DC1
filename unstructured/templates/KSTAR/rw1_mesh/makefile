all : input
	m3dc1_mfmgen input

fetch :
	wget -O kstar_rw1.txt 'https://drive.google.com/uc?export=download&id=1yDrezaujSYZuuQ9nwjKEMzkjDI-Eqjgn'
	wget -O kstar_rw1-7K0.smb 'https://drive.google.com/uc?export=download&id=1ngNJNHIroz1pKk5ntQYsDgjZzUNj7Bci'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
