all : input
	m3dc1_mfmgen input

fetch :
	wget -O nstx_fw1.txt 'https://drive.google.com/uc?export=download&id=1ON0LEcIdZ_6XzuDRebL4tfcl-6pb3NLR'
	wget -O nstx_fw1-14K0.smb 'https://drive.google.com/uc?export=download&id=1eXd4bEPJDID9M32hOnnN-_9QDXWo1gvb'

.PHONY : clean all fetch

clean :
	rm -fr *.smb *.sms *.smd *.dmg *.vtk *.txt *.log plot_geo *~
