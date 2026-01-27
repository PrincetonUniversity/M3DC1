all : input
	m3dc1_mfmgen input

fetch :
	wget -O aug_rw1.txt 'https://drive.google.com/uc?export=download&id=1f7BQJpc8rtrQloeCuiWF7nJkODZK-ND2'
	wget -O aug_rw1-7K0.smb 'https://drive.google.com/uc?export=download&id=1ByuZbgTikomcyKgYZUNVzb9dU6QUlgg1'

.PHONY : clean all fetch

clean :
	rm -fr *~ *.smb *.smd *.sms *.dmg *.vtk *.txt *.log plot_geo
