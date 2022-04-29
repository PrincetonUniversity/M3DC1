#!/usr/bin/env python3
#
# injection_rate: Calculate rate of injected impurities
#
# Coded on 03/23/2022 by:
# Andreas Kleiner:    akleiner@pppl.gov

import numpy as np
import math
import fpy
from m3dc1.read_h5 import readParameter
from m3dc1.read_h5 import readC1File

def injection_rate(sim=None,filename='C1.h5'):
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    h5file = sim._all_attrs
    
    time_n,kprad_n = readC1File(scalar='kprad_n')
    time_n0,kprad_n0 = readC1File(scalar='kprad_n0')
    
    B0 = readParameter('b0_norm',h5file=h5file)
    N0 = readParameter('n0_norm',h5file=h5file)
    L0 = readParameter('l0_norm',h5file=h5file)
    mi = readParameter('ion_mass',h5file=h5file)
    mu0=1.0 # Magnetic permeability in cgs units
    mp = 1.67E-24
    
    V0 = B0/(4*math.pi*mu0*mi*mp*N0)**0.5
    T0 = L0/V0
    
    #print(N0*1E6)
    #print(T0)
    #print(time_n[1]-time_n[0])
    #print(kprad_n[1]-kprad_n[0])
    rate_n = ((kprad_n[1]-kprad_n[0])*N0*1E6)/((time_n[1]-time_n[0])*T0)
    rate_n0 = ((kprad_n0[1]-kprad_n0[0])*N0*1E6)/((time_n0[1]-time_n0[0])*T0)
    
    print('Rate of total impurities: '+str(rate_n))
    print('Rate of neutral impurities: '+str(rate_n0))
    return
