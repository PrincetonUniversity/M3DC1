

.PHONY : clean interpolate

interpolate :
	cp g158115.04701_590 geqdsk	
	cp ${M3DC1_DIR}/device_data/DIII-D/coil.dat .
	a2cc a158115.04701_590 > current.dat
	extract_profiles.sh p158115.04701_590
	cp ../Create_Multiregion_Mesh/diiid-8K0.smb .
	cp ../Create_Multiregion_Mesh/diiid.txt .
	sbatch batchjob.$M3DC1_ARCH.interpolate

clean :
	rm -fr *~ slurm* *.h5 part*.smb ts0-adapted*.smb geqdsk coil.dat current.dat* mesh0.smb offset.txt debug-out eqdsk-out testout profile* C1ke interpolated normcurv
