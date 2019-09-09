# Converts from M3DC1 units to mks and vice versa
# Makes use of the pint module
#
#

import pint
import os

def unit_conv(array, arr_dim='M3DC1', time=0, length=0, particles=0, magnetic_field=0, current=0, current_density=0, diffusion=0, energy=0, force=0, pressure=0, resistivity=0, temperature=0, velocity=0, voltage=0, viscosity=0, thermal_conductivity=0, electric_field=0):
    """
    Converts an array from M3DC1 units to mks or vice versa. arr_dim
    contains the type of dimension the array is in (so 'M3DC1', or 
    'mks'). The routine will return a new array in the opposing 
    dimensions. One must specify in what dimensions the array is given. 
    If, for example, an array is given in units of energy**2/current_density,
    input energy=2 and current_density=-1.
    """
    
    dim_path = os.environ['DIM_TXT']
    ureg = pint.UnitRegistry(dim_path)
    
    if arr_dim=='M3DC1':
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
                        ureg.Pascal**pressure * \
                        ureg.OhmMeter**resistivity * \
                        ureg.eV**temperature * \
                        ureg.MeterPerSecond**velocity * \
                        ureg.Volts**voltage * \
                        ureg.KilogramPerMeterPerSecond**viscosity * \
                        ureg.PerMeterPerSecond**thermal_conductivity * \
                        ureg.VoltsPerMeter**electric_field )
                        
    if arr_dim=='mks':
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
                        ureg.M3DC1pressure**pressure * \
                        ureg.M3DC1resistivity**resistivity * \
                        ureg.M3DC1temperature**temperature * \
                        ureg.M3DC1velocity**velocity * \
                        ureg.M3DC1voltage**voltage * \
                        ureg.M3DC1viscosity**viscosity * \
                        ureg.M3DC1thermalconductivity**thermal_conductivity * \
                        ureg.M3DC1electricfield**electric_field )

    return array_conv.magnitude
