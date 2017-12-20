#!/bin/bash

echo "!simVer: 0   - the latest SimModeler (5.0 as of 06/15/16)" > input
echo "!        100 - SimModeler 5.0" >> input
echo "!        90  - SimModeler 4.0" >> input
echo "simVer 90" >> input
echo "!model type: 0 - parameterized vacuum" >> input
echo "!            1 - piece-wise linear vacuum; " >> input
echo "!            2 - piece-wise polynomial vacuum" >> input
echo "!            3 - three-regions" >> input
echo "modelType 0" >> input
echo "outFile analytic" >> input
echo "meshSize $2" >> input
echo "" >> input
echo "! for parametrized vacuum model only." >> input
echo "useVacuumParams 1" >> input
echo "vacuumParams" `head -n 1 $1` >> input
echo "numVacuumPts 20" >> input
echo "" >> input
echo "! 1 to multiply coord and params of nodes on vacuum boundary by vacuumFactor" >> input
echo "! If vacuumFactor is not specified, the default is 2*PI" >> input
echo "! If adjustVacuumParams=1, using SimModeler for further meshing is not supported and simmmetrix mesh file is not generated" >> input
echo "adjustVacuumParams 0" >> input
echo "vacuumFactor  6.28319" >> input

m3dc1_meshgen input
