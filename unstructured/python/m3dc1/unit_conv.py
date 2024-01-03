#!/usr/bin/env python3
# Converts from M3DC1 units to mks and vice versa
# Makes use of the pint module
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import pint
import os
import math
import fpy
from m3dc1.read_h5 import readParameter

def unit_conv(array, arr_dim='M3DC1', filename='C1.h5', sim=None, time=0, length=0, particles=0, magnetic_field=0, current=0, current_density=0, diffusion=0, energy=0, force=0, magnetic_flux=0, pressure=0, resistivity=0, temperature=0, velocity=0, voltage=0, viscosity=0, thermal_conductivity=0, electric_field=0):
    """
    Converts an array from M3DC1 units to mks or vice versa. arr_dim
    contains the type of dimension the array is in (so 'M3DC1', or 
    'mks'). The routine will return a new array in the opposing 
    dimensions. One must specify in what dimensions the array is given. 
    If, for example, an array is given in units of energy**2/current_density,
    input energy=2 and current_density=-1.
    """
    
    ureg = pint.UnitRegistry(None)
    ureg.define('second       = [time] = s = sec')
    ureg.define('meter 	     = [length]')
    ureg.define('particles      = [number_of_particles] = N')
    ureg.define('Gauss = (1.0*10.0**-4) * Tesla')
    ureg.define('Tesla = [magnetic_field_strength] = T')
    ureg.define('Ampere       = [current] = A = Amp')
    ureg.define('AmperePerSquareMeter = [current_density]')
    ureg.define('SquareMeterPerSecond = [diffusion]')
    ureg.define('Joules = [energy] = J')
    ureg.define('Newton = [force] = N')
    ureg.define('Weber = [magnetic_flux] = Wb')
    ureg.define('Pascal = [pressure] = Pa')
    ureg.define('OhmMeter = [resistivity]')
    ureg.define('eV = [temperature]')
    ureg.define('MeterPerSecond = [velocity]')
    ureg.define('Volts        = [voltage]')
    ureg.define('KilogramPerMeterPerSecond = [viscosity]')
    ureg.define('PerMeterPerSecond 	 = [thermal_conductivity]')
    ureg.define('VoltsPerMeter      = [electric_field]')

    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    h5file = sim._all_attrs

    B0 = readParameter('b0_norm',h5file=h5file)
    N0 = readParameter('n0_norm',h5file=h5file)
    L0 = readParameter('l0_norm',h5file=h5file)
    mi = readParameter('ion_mass',h5file=h5file)
    version = readParameter('version',h5file=h5file)
    if version >= 23:
        Zeff = readParameter('z_ion',h5file=h5file)
    else:
        Zeff = readParameter('zeff',h5file=h5file)
    
    #mi=1 # Uncomment for benchmarking
    
    mu0=1.0 # Magnetic permeability in cgs units
    mp = 1.67E-24
    c = 3.00E10 # Vacuum speed of light
    e = 1.602176634E-12 # Elementary charge
    # Convert to mks units, remove below stuff eventually
    #mu0 = 4.0*math.pi*1.0E-7
    #B0=B0/10000.0
    #N0=N0*1.0E6
    #L0=L0/100.0

    V0 = B0/(4*math.pi*mu0*mi*mp*N0)**0.5
    T0 = L0/V0
    I0 = c*B0*L0/(4.0*math.pi) #For current normalization
    J0 = c*B0/(4.0*math.pi*L0)
    diffus0 = L0**2/T0
    E0 = V0*B0/c
    energy0 = B0**2*L0**3 / (4.0*math.pi)
    force0 = B0**2*L0**2 / (4.0*math.pi)
    magflux0 = B0*L0**2
    pressure0 = B0**2 / (4.0*math.pi)
    eta0 = 4.0*math.pi*T0*V0**2 / c**2
    temp0 = B0**2/(4.0*math.pi*N0*e)
    thcond0 = N0*L0**2/T0
    visc0 = mi*mp*N0*L0**2/T0
    volt0 = L0*V0*B0/c
    len0 = L0
    
    #Now everything is in cgs units, next convert to mks
    V0 = V0/1.0E2
    I0 = I0/3.0E9
    J0 = J0/3.0E5
    diffus0 = diffus0/1.0E4
    E0 = E0/(1/3*1E-4)
    energy0 = energy0/1.0E7
    force0 = force0/1.0E5
    magflux0 = magflux0/1.0E8
    pressure0 = pressure0/10
    eta0 = eta0/(1/9*1.0E-9)
    temp0 = temp0
    thcond0 = thcond0*100
    visc0 = visc0/10
    volt0 = volt0/(1/3*1.0E-2)
    len0 = len0/100.0

    
    #print('V0: '+'{:.5E}'.format(V0))
    #print('T0: '+'{:.5E}'.format(T0))
    #print('I0: '+'{:.5E}'.format(I0))
    #print('J0: '+'{:.5E}'.format(J0))
    #print('diffus0: '+'{:.5E}'.format(diffus0))
    #print('E0: '+'{:.5E}'.format(E0))
    #print('energy0: '+'{:.5E}'.format(energy0))
    #print('force0: '+'{:.5E}'.format(force0))
    #print('pressure0: '+'{:.5E}'.format(pressure0))
    #print('eta0: '+'{:.5E}'.format(eta0))
    #print('temp0: '+'{:.5E}'.format(temp0))
    #print('thcond0: '+'{:.5E}'.format(thcond0))
    #print('visc0: '+'{:.5E}'.format(visc0))
    #print('volt0: '+'{:.5E}'.format(volt0))
    
    ureg.define('M3DC1time    = ('+str(T0)+') * second')
    ureg.define('M3DC1velocity = ('+str(V0)+') * MeterPerSecond')
    ureg.define('M3DC1length  = '+str(len0)+' * meter')
    ureg.define('M3DC1particles = (1.00*10.0**20) * particles')
    ureg.define('M3DC1magneticfield = (1.00*10.0**4) * Gauss')
    ureg.define('M3DC1current = ('+str(I0)+') * Ampere')
    ureg.define('M3DC1currentdensity = ('+str(J0)+') * AmperePerSquareMeter')
    ureg.define('M3DC1diffusion = ('+str(diffus0)+') * SquareMeterPerSecond')
    ureg.define('M3DC1energy = ('+str(energy0)+') * Joules')
    ureg.define('M3DC1force = ('+str(force0)+') * Newton')
    ureg.define('M3DC1magflux = ('+str(magflux0)+') * Weber')
    ureg.define('M3DC1pressure = ('+str(pressure0)+') * Pascal')
    ureg.define('M3DC1resistivity = ('+str(eta0)+') * OhmMeter')
    ureg.define('M3DC1temperature = ('+str(temp0)+') * eV')
    ureg.define('M3DC1velocity = ('+str(V0)+') * MeterPerSecond')
    ureg.define('M3DC1voltage = ('+str(volt0)+') * Volts')
    ureg.define('M3DC1viscosity = ('+str(visc0)+') * KilogramPerMeterPerSecond')
    ureg.define('M3DC1thermalconductivity = ('+str(thcond0)+') * PerMeterPerSecond')
    ureg.define('M3DC1electricfield = ('+str(E0)+') * VoltsPerMeter')
    
    if arr_dim.lower()=='m3dc1':
        array_dimfull = array * \
                        ureg.M3DC1time**time * \
                        ureg.M3DC1length**length * \
                        ureg.M3DC1particles**particles * \
                        ureg.M3DC1magneticfield**magnetic_field * \
                        ureg.M3DC1current**current * \
                        ureg.M3DC1currentdensity**current_density * \
                        ureg.M3DC1diffusion**diffusion * \
                        ureg.M3DC1energy**energy * \
                        ureg.M3DC1force**force * \
                        ureg.M3DC1magflux**magnetic_flux * \
                        ureg.M3DC1pressure**pressure * \
                        ureg.M3DC1resistivity**resistivity * \
                        ureg.M3DC1temperature**temperature * \
                        ureg.M3DC1velocity**velocity * \
                        ureg.M3DC1voltage**voltage * \
                        ureg.M3DC1viscosity**viscosity * \
                        ureg.M3DC1thermalconductivity**thermal_conductivity * \
                        ureg.M3DC1electricfield**electric_field
        array_conv   =  array_dimfull.to( \
                        ureg.sec**time * \
                        ureg.meter**length * \
                        ureg.particles**particles * \
                        ureg.Tesla**magnetic_field * \
                        ureg.Ampere**current * \
                        ureg.AmperePerSquareMeter**current_density * \
                        ureg.SquareMeterPerSecond**diffusion * \
                        ureg.Joules**energy * \
                        ureg.Newton**force * \
                        ureg.Weber**magnetic_flux * \
                        ureg.Pascal**pressure * \
                        ureg.OhmMeter**resistivity * \
                        ureg.eV**temperature * \
                        ureg.MeterPerSecond**velocity * \
                        ureg.Volts**voltage * \
                        ureg.KilogramPerMeterPerSecond**viscosity * \
                        ureg.PerMeterPerSecond**thermal_conductivity * \
                        ureg.VoltsPerMeter**electric_field )
                        
    if arr_dim.lower()=='mks':
        array_dimfull = array * \
                        ureg.sec**time * \
                        ureg.meter**length * \
                        ureg.particles**particles * \
                        ureg.Tesla**magnetic_field * \
                        ureg.Ampere**current * \
                        ureg.AmperePerSquareMeter**current_density * \
                        ureg.SquareMeterPerSecond**diffusion * \
                        ureg.Joules**energy * \
                        ureg.Newton**force * \
                        ureg.Weber**magnetic_flux * \
                        ureg.Pascal**pressure * \
                        ureg.OhmMeter**resistivity * \
                        ureg.eV**temperature * \
                        ureg.MeterPerSecond**velocity * \
                        ureg.Volts**voltage * \
                        ureg.KilogramPerMeterPerSecond**viscosity * \
                        ureg.PerMeterPerSecond**thermal_conductivity * \
                        ureg.VoltsPerMeter**electric_field
        array_conv   =  array_dimfull.to( \
                        ureg.M3DC1time**time * \
                        ureg.M3DC1length**length * \
                        ureg.M3DC1particles**particles * \
                        ureg.M3DC1magneticfield**magnetic_field * \
                        ureg.M3DC1current**current * \
                        ureg.M3DC1currentdensity**current_density * \
                        ureg.M3DC1diffusion**diffusion * \
                        ureg.M3DC1energy**energy * \
                        ureg.M3DC1force**force * \
                        ureg.M3DC1magflux**magnetic_flux * \
                        ureg.M3DC1pressure**pressure * \
                        ureg.M3DC1resistivity**resistivity * \
                        ureg.M3DC1temperature**temperature * \
                        ureg.M3DC1velocity**velocity * \
                        ureg.M3DC1voltage**voltage * \
                        ureg.M3DC1viscosity**viscosity * \
                        ureg.M3DC1thermalconductivity**thermal_conductivity * \
                        ureg.M3DC1electricfield**electric_field )

    return array_conv.magnitude
