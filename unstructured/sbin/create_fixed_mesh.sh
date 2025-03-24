#!/bin/bash -l

meshgen=$(which m3dc1_meshgen)

if [ -x "${meshgen}" ]; then
    echo "m3dc1_meshgen found."
    FILE="input_AnalyticModel"
    rm -f $FILE

    echo "!simVer: 0   - the latest SimModeler (5.0 as of 06/15/16)" > $FILE
    echo "!        100 - SimModeler 5.0" >> $FILE
    echo "!        90  - SimModeler 4.0" >> $FILE
    echo "simVer 90" >> $FILE
    echo "!model type: 0 - parameterized vacuum" >> $FILE
    echo "!            1 - piece-wise linear vacuum; " >> $FILE
    echo "!            2 - piece-wise polynomial vacuum" >> $FILE
    echo "!            3 - three-regions" >> $FILE
    echo "modelType 0" >> $FILE
    echo "outFile analytic" >> $FILE
    echo "meshSize $2" >> $FILE
    echo "" >> $FILE
    echo "! for parametrized vacuum model only." >> $FILE
    echo "useVacuumParams 1" >> $FILE
    echo "vacuumParams" `head -n 1 $1` >> $FILE
    echo "numVacuumPts 20" >> $FILE
    echo "" >> $FILE
    echo "! 1 to multiply coord and params of nodes on vacuum boundary by vacuumFactor" >> $FILE
    echo "! If vacuumFactor is not specified, the default is 2*PI" >> $FILE
    echo "! If adjustVacuumParams=1, using SimModeler for further meshing is not supported and simmmetrix mesh file is not generated" >> $FILE
    echo "adjustVacuumParams 0" >> $FILE
    echo "vacuumFactor  6.28319" >> $FILE

    m3dc1_meshgen $FILE
else
    echo "m3dc1_meshgen not found.  Trying create_smb"

    if [ ! -e seed0.smb ]; then
	cp $M3DC1_DIR/bin/seed0.smb .
    fi

    srun -n 1 create_smb $1 $2
fi
