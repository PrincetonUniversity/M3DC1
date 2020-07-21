#!/bin/bash

B0=1e4
N0=1e14
L0=100.
MU=1
ZEFF=1
PI=3.14159265359
C=2.9979e10
MP=1.6726e-24
E=4.8032e-10

# Check if C1input exists
if [ ! -f C1input ]; then
    echo 'WARNING: C1input file not found.  Using default values.'
else

    # Parse C1input for normalizations
    BLINE=`grep b0_norm C1input`
    if [ $? == 0 ]; then
	B0=`echo $BLINE | awk '{ print $3 }'`
    else
	echo "b0_norm not found.  Using default value"
    fi
    NLINE=`grep n0_norm C1input`
    if [ $? == 0 ]; then
	N0=`echo $NLINE | awk '{ print $3 }'`
    else
	echo "n0_norm not found.  Using default value"
    fi
    LLINE=`grep l0_norm C1input`
    if [ $? == 0 ]; then
	L0=`echo $LLINE | awk '{ print $3 }'`
    else
	echo "l0_norm not found.  Using default value"
    fi
    MULINE=`grep ion_mass C1input`
    if [ $? == 0 ]; then
	MU=`echo $MULINE | awk '{ print $3 }'`
    else
	echo "ion_mass not found.  Using default value"
    fi
    ZLINE=`grep 'zeff\|z_ion' C1input`
    if [ $? == 0 ]; then
	ZEFF=`echo $ZLINE | awk '{ print $3 }'`
    else
	echo "zeff not found.  Using default value"
    fi
fi

B0_MKS=$(echo | awk "{ print $B0/1.e4 }" )
N0_MKS=$(echo | awk "{ print $N0*1.e6 }" | awk '{ printf "%.5e", $1 }' )
L0_MKS=$(echo | awk "{ print $L0/100. }" )

V0=$(echo - | awk "{ print $B0/sqrt(4.*$PI*$MU*$MP*$N0) }" )
V0_MKS=$(echo | awk "{ print $V0/100. }" )

P0=$(echo - | awk "{ print $B0*$B0/(4.*$PI) }" )
P0_MKS=$(echo - | awk "{ print $P0/10. }" )

T0=$(echo - | awk "{ print $L0/$V0 }" )
FREQ0=$(echo - | awk "{ print 1./$T0 }" )
AFREQ0=$(echo - | awk "{ print 2*$PI*$FREQ0 }" )

W0=$(echo - | awk "{ print $P0/$N0 }" )
W0_MKS=$(echo - | awk "{ print $W0/1e7 }" )
W0_EV=$(echo - | awk "{ print $W0/1.6022e-12 }" )

J0=$(echo - | awk "{ print $C*$B0/(4.*$PI*$L0) }" )
J0_MKS=$(echo - | awk "{ print $J0 / 3e5 }" )

ETA0=$(echo - | awk "{ print 4*$PI*$T0*$V0*$V0/($C*$C) }" )
ETA0_MKS=$(echo - | awk "{ print $ETA0*$C*$C/1e11 }" )

MU0=$(echo - | awk "{ print $MP*$MU*$N0*$L0*$L0/$T0 }" )
MU0_MKS=$(echo - | awk "{ print $MU0/10. }" )

KAP0=$(echo - | awk "{ print $N0*$L0*$L0/$T0 }" | awk '{ printf "%.5e", $1 }' )
KAP0_MKS=$(echo - | awk "{ print 100.*$N0*$L0*$L0/$T0 }" | awk '{ printf "%.5e", $1 }' )

VOLT0=$(echo - | awk "{ print $V0*$B0*$L0/$C }" )
VOLT0_MKS=$(echo - | awk "{ print $V0*$B0*$L0/1e8 }" | awk '{ printf "%.5e", $1 }')


echo "ion mass = $MU m_p"
echo "Zeff = $ZEFF"
echo "B0 = $B0 G = $B0_MKS (T)"
echo "n0 = $N0 /cm^3 = $N0_MKS particles/m^3"
echo "L0 = $L0 cm = $L0_MKS  m"
echo "t0 = $T0 s"
echo "1/t0 =           $FREQ0 Hz"
echo "2pi/t0           $AFREQ0 rad/s"
echo "Velocity:        $V0 cm/s = $V0_MKS m/s"
echo "Pressure:        $P0 dyne/cm^2 = $P0_MKS pascal"
echo "Temperature:     $W0_EV eV"
echo "Current density: $J0 s.a./cm^2 = $J0_MKS A/m^2"
echo "Resistivity:     $ETA0 s = $ETA0_MKS Ohm-m"
echo "Viscosity:       $MU0 g/(cm-s) = $MU0_MKS kg/(m-s)"
echo "Thermal diff.:   $KAP0 particles/(cm-s) = $KAP0_MKS particles/(m-s)"
echo "Voltage:         $VOLT0 s.v. = $VOLT0_MKS V"
