#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: akleiner
"""
import fpy
import numpy as np
import math
import glob
import os
from termcolor import colored
import matplotlib.pyplot as plt
from m3dc1.unit_conv  import unit_conv

#-------------------------------------------
# Plotting and formatting routines
#-------------------------------------------

def plot2d(xdata,ydata,equal=False):
    """
    Create a simple 2D line plot.

    Arguments:

    **xdata**
    Array containing x values

    **ydata**
    Array containing y values
    
    **equal**
    If True, x and y axes have same scaling
    """
    plt.figure()
    plt.plot(xdata,ydata)
    plt.grid(True)
    plt.show()
    if equal:
        plt.axis('equal')
    
    return


# Determines reasonables values for lower and upper axis limits.
def get_axlim(val,minmax,delta=0.1):
    if minmax=='max':
        lim = math.ceil(val*2.)/2.
        if np.abs(lim - val)<delta:
            lim = lim + delta
    elif minmax=='min':
        lim = math.floor(val*2.)/2.
        if np.abs(lim - val)<delta:
            lim = lim - delta
    return lim



# Formats float to be in engineering notation
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \cdot 10^{{{}}}$'.format(a, b)


#-------------------------------------------
# Read simulations via fusion-io
#-------------------------------------------
# sets up arrays for sim and time
def setup_sims(sim,filename,time,linear,diff):
    # make iterable
    if not isinstance(sim,fpy.sim_data):
        sim = np.empty(0)
        filename = np.atleast_1d(filename)
        for f in filename:
            sim = np.append(sim,fpy.sim_data(f))
    else:
        sim = np.atleast_1d(sim)

    time = np.atleast_1d(time)

    if len(time)==1 and len(sim)>1:
        time = np.repeat(time,len(sim))
    elif len(sim)==1 and len(time)>1:
        sim = np.repeat(sim,len(time))
    elif len(sim) != len(time):
        raise RuntimeError('Length of time does not match length of sim/filename')

    if linear:
        if len(sim)>1:
            raise RuntimeError('Provide a single simulation for linear=True')
        if (time[0]==-1) or ((time[0] is None) and (sim[0].timeslice==-1)):
            raise RuntimeError('time or sim.timeslice must be greater than -1 for linear=True')
        sim = np.repeat(sim,2)
        time = np.append(time,-1)

    ### Input error handling ###
    if diff and (len(time) != 2):
        raise RuntimeError('Please input twoxs times for differences or specify two sim_data objects.')

    if diff and linear:
        raise RuntimeError('Please choose diff or linear, not both.')

    if (not diff) and (not linear) and (len(sim)>1):
        raise RuntimeError('Multiple simulations detected. Please set diff=True or input single slices')

    return sim, time

# identify index for given coordinate
def get_field_idx(coord):

    field_idx = {'R':0, 'scalar':0, 'phi':1, 'Z':2}
    if coord in field_idx:
        return field_idx[coord]
    elif coord in ['poloidal', 'radial', 'vector', 'tensor']:
        return None
    else:
        raise RuntimeError('Please enter valid coordinate. Accepted: \'R\', \'phi\', \'Z\', \'poloidal\', \'radial\', \'scalar\', \'vector\'')





def read_floats(string,length=16):
    """
    Converts string containing multiple numbers in scientific notation to a list of floats
    
    Arguments:

    **string**
    String to convert
    
    **length**
    Number of characters corresponding to one floating point number,
    e.g. -0.12345E-01 has length=12
    """
    float_list = [string[start:start+length] for start in range(0, len(string), length)]
    if float_list[-1] == '\n':
        float_list = float_list[:-1]
    for i in range(len(float_list)):
        float_list[i] = float(float_list[i])
    return float_list



#-------------------------------------------
# Mathematical and numerical routines
#-------------------------------------------

def deriv(y,x=None):
    """
    Calculates derivative of an array y using three-point (quadratic) Lagrangian interpolation.
    This function resembles the corresponding IDL function deriv(). Note that the arguments
    x and y are switched.

    Arguments:

    **y**
    Array containing the function values y

    **x**
    Array containing x locations where y is given. If x is omitted, y is assumed to
    be evenly spaced
    """
    yp = np.zeros_like(y)
    
    if isinstance(x, (np.ndarray,list)):
        if len(x)!=len(y):
            raise Exception('x and y do not have the same length.')
        for i in range(len(y)):
            if i==0:
                x01 = x[0]-x[1]
                x02 = x[0]-x[2]
                x12 = x[1]-x[2]
                yp[i] = y[0]*(x01 + x02)/(x01*x02) - y[1]*x02 / (x01*x12) + y[2]*x01 / (x02*x12)
            elif i == len(y)-1:
                x01 = x[-3]-x[-2]
                x02 = x[-3]-x[-1]
                x12 = x[-2]-x[-1]
                yp[i] = -y[-3]*x12/(x01*x02) + y[-2]*x02/(x01*x12) - y[-1]*(x02 + x12)/(x02*x12)
            else:
                x01 = x[i-1]-x[i]
                x02 = x[i-1]-x[i+1]
                x12 = x[i]-x[i+1]
                yp[i] = y[i-1]*x12/(x01*x02) + y[i]*(1.0/x12 - 1.0/x01) - y[i+1]*x01/(x02*x12)
    else:
        #Derivative for evenly-spaced array, used when x is not provided:
        for i in range(len(y)):
            if i==0:
                yp[i] = (-3.0*y[0] + 4.0*y[1] - y[2])/2.0
            elif i == len(y)-1:
                yp[i] = (3.0*y[-1] - 4.0*y[-2] + y[-3])/2.0
            else:
                yp[i] = (y[i+1] - y[i-1])/2.0
    return yp



def smooth(vin, w, nan='replace'): 
    """
    Calculates the boxcar average of an array. Closely resembles the IDL smooth function.

    Arguments:

    **vin**
    Input array, can be 1D or 2D

    **w**
    width of smoothing window

    **nan**
    Choose how to treat NaNs.
    'replace': Ignore NaNs and consider only finite values
    'propagate': Keep NaNs
    """
    # create output array:
    vout=np.copy(vin)

    # If w is even, add 1
    if w % 2 == 0:
        w = w + 1
        print('Window width is even. Setting w=w+1='+str(w))

    # get the size of each dim of the input:
    dims = len(vin.shape)
    # Check for dimension af array
    if dims==1:
        r = vin.shape[0]
        c = 0
    elif dims==2:
        r, c = vin.shape
    else:
        raise Exception('Please provide either a one-dimensional or two-dimensional array!')
    
    # Assume that the width of the window w is always square.
    startrc = int((w - 1)/2)
    stopr = int(r - ((w + 1)/2) + 1)
    stopc = int(c - ((w + 1)/2) + 1)

    if dims==1:
        for row in range(startrc,stopr):
            # Determine the window
            startwr = int(row - (w/2))
            stopwr = int(row + (w/2) + 1)
            window = vin[startwr:stopwr]
            if nan == 'replace':
                # If we're replacing Nans, then select only the finite elements
                window = window[np.isfinite(window)]
            # Calculate the mean of the window
            vout[row] = np.mean(window)
    elif dims == 2:
        for col in range(startrc,stopc):
            # Start and stop indices for column
            startwc = int(col - (w/2))
            stopwc = int(col + (w/2) + 1)
            for row in range (startrc,stopr):
                # Determine the window
                startwr = int(row - (w/2))
                stopwr = int(row + (w/2) + 1)
                window = vin[startwr:stopwr, startwc:stopwc]
                if nan == 'replace':
                    # Ignore NaNs and only consider finite values
                    window = window[np.isfinite(window)]
                # Calculate the mean inside the window
                vout[row,col] = np.mean(window)
    return vout



def PolygonArea(x,y):
    """
    Calculates the area inside a polygon

    Arguments:

    **x**
    Array of polygon point x values

    **y**
    Array of polygon point x values
    """
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

#-------------------------------------------
# Unit conversion
#-------------------------------------------



def get_unitexpns():
    return {'time':0, 'length':0, 'particles':0, 'magnetic_field':0,
            'current':0, 'current_density':0, 'diffusion':0, 'energy':0,
            'force':0, 'pressure':0, 'resistivity':0, 'temperature':0,
            'velocity':0, 'voltage':0, 'viscosity':0,
            'thermal_conductivity':0, 'electric_field':0}


# Returns field label depending on chosen system of units
def get_fieldlabel(units,field,shortlbl=False):

    labels = {'j':'current density', 'ni':'ion density','ne':'electron density',
              'v':'velocity', 'B':'magnetic field strength', 'p':'pressure',
              'pi':'ion pressure', 'pe':'electron pressure',
              'ti':'ion temperature', 'te':'electron temperature',
              'A':'vector potential', 'gradA':'grad vector potential',
              'E':'electric field', 'default':field}
    
    short_labels = {'j':'j', 'ni':'$n_{i}$','ne':'$n_{e}$',
              'v':'v', 'B':'B', 'p':'p',
              'pi':'$p_{i}$', 'pe':'$p_{e}$',
              'ti':'$T_{i}$', 'te':'$T_{e}$',
              'A':'A', 'gradA':'$grad A$ ',
              'E':'E', 'default':field}

    if units.lower()=='m3dc1':
        units = {'default':'M3DC1 units'}
    elif units.lower()=='mks':
        units = {'j':'$A/m^2$', 'ni':'particles/$m^3$', 'ne':'particles/$m^3$',
                 'v':'$m/s$', 'B':'$T$', 'p':'$Pa$', 'pi':'$Pa$', 'pe':'$Pa$',
                 'ti':'$eV$', 'te':'$eV$', 'A':'$Tesla \cdot m$',
                 'gradA':'$(Tesla \cdot m)$ / (m or rad)', 'E':'$V/m$',
                 'default':'MKS units'}

    if field in labels:
        label = short_labels[field] if shortlbl else labels[field]
    else:
        label = labels['default']

    if field in units:
        unit = units[field]
    else:
        unit = units['default']

    return label, unit


def get_conv_field(units,field,field1_ave,filename='C1.h5',sim=None):
    """
    Returns converted field depending on chosen system of units
    """

    expns = get_unitexpns()
    fields = {'j':{'current_density':1}, 'ni':{'particles':1,'length':-3},
              'ne':{'particles':1,'length':-3}, 'v':{'velocity':1},
              'B':{'magnetic_field':1}, 'p':{'pressure':1}, 'pi':{'pressure':1},
              'pe':{'pressure':1}, 'ti':{'temperature':1},
              'te':{'temperature':1}, 'A':{'magnetic_field':1,'length':1},
              'E':{'electric_field':1},
              }
    if field in fields:
        expns.update(fields[field])

    if units.lower()=='m3dc1':
        if not isinstance(sim,fpy.sim_data):
            sim = fpy.sim_data(filename=filename)
        field1_ave = unit_conv(field1_ave,arr_dim='mks',sim=sim,**expns)
    return field1_ave


def get_conv_trace(units,trace,trace_arr,filename='C1.h5',sim=None,itor=1,custom=None):
    """
    Returns converted time trace depending on chosen system of units
    """

    expns = get_unitexpns()

    traces = {'Ave_P':{'pressure':1}, 'E_K3':{'energy':1},
              'E_K3D':{'energy':1,'time':-1}, 'E_K3H':{'energy':1,'time':-1},
              'E_KP':{'energy':1}, 'E_KPD':{'energy':1,'time':-1},
              'E_KPH':{'energy':1,'time':-1}, 'E_KT':{'energy':1},
              'E_KTD':{'energy':1,'time':-1}, 'E_KTH':{'energy':1,'time':-1},
              'E_MP':{'energy':1}, 'E_MPD':{'energy':1,'time':-1},
              'E_MPH':{'energy':1,'time':-1}, 'E_MT':{'energy':1},
              'E_MTD':{'energy':1,'time':-1}, 'E_MTH':{'energy':1,'time':-1},
              'E_P':{'energy':1}, 'E_PD':{'energy':1,'time':-1},
              'E_PE':{'energy':1}, 'E_PH':{'energy':1,'time':-1},
              'E_grav':{'energy':1}, 'Flux_kinetic':{'energy':1,'time':-1},
              'Flux_poynting':{'energy':1,'time':-1},
              'Flux_pressure':{'energy':1,'time':-1},
              'Flux_thermal':{'energy':1,'time':-1}, 'IP_co':{'current':1},
              'IP_sn':{'current':1}, 'M_IZ':{'current':1,'length':2},
              'M_IZ_co':{'current':1,'length':2},
              'M_IZ_sn':{'current':1,'length':2},
              'Parallel_viscous_heating':{'energy':1,'time':-1},
              'Particle_Flux_convective':{'particles':1,'time':-1},
              'Particle_Flux_diffusive':{'particles':1,'time':-1},
              'Particle_source':{'particles':1,'time':-1},
              'Torque_com':{'force':1,'length':1},
              'Torque_em':{'force':1,'length':1},
              'Torque_gyro':{'force':1,'length':1},
              'Torque_parvisc':{'force':1,'length':1},
              'Torque_sol':{'force':1,'length':1},
              'Torque_visc':{'force':1,'length':1},
              'W_M':{'energy':1}, 'W_P':{'energy':1},
              'Wall_Force_n0_x':{'force':1}, 'Wall_Force_n0_x_halo':{'force':1},
              'Wall_Force_n0_y':{'force':1}, 'Wall_Force_n0_z':{'force':1},
              'Wall_Force_n0_z_halo':{'force':1}, 'Wall_Force_n1_x':{'force':1},
              'Wall_Force_n1_y':{'force':1},
              'angular_momentum':{'force':1,'length':1,'time':1},
              'angular_momentum_p':{'force':1,'length':1,'time':1},
              'area':{'length':2}, 'area_p':{'length':2},
              'brem_rad':{'energy':1,'time':-1},
              'circulation':{'velocity':1,'length':1}, 'dt':{'time':1},
              'electron_number':{'particles':1},
              'helicity':{'magnetic_field':2,'length':4},
              'i_control%err_i':{'current':1,'time':1},
              'i_control%err_p_old':{'current':1},
              'ion_loss':{'energy':1,'time':-1},
              'line_rad':{'energy':1,'time':-1},
              'loop_voltage':{'voltage':1},
              'n_control%err_i':{'particles':1,'time':1},
              'n_control%err_p_old':{'particles':1},
              'particle_number':{'particles':1},
              'particle_number_p':{'particles':1},
              'psi0':{'magnetic_field':1,'length':2},
              'psi_lcfs':{'magnetic_field':1,'length':2},
              'psimin':{'magnetic_field':1,'length':2},
              'radiation':{'energy':1,'time':-1},
              'reconnected_flux':{'magnetic_field':1,'length':1+itor},
              'reck_rad':{'energy':1,'time':-1},
              'recp_rad':{'energy':1,'time':-1}, 'runaways':{'particles':1},
              'temax':{'temperature':1}, 'time':{'time':1},
              'toroidal_current':{'current':1},
              'toroidal_current_p':{'current':1},
              'toroidal_current_w':{'current':1},
              'toroidal_flux':{'magnetic_field':1,'length':2},
              'toroidal_flux_p':{'magnetic_field':1,'length':2},
              'volume':{'length':3}, 'volume_p':{'length':3}, 'xmag':{'length':1},
              'xnull':{'length':1}, 'xnull2':{'length':1}, 'zmag':{'length':1},
              'znull':{'length':1}, 'znull2':{'length':1},
              'bharmonics':{'energy':1}, 'keharmonics':{'energy':1},
              'cauchy_fraction':{}, 'cloud_pel':{}, 'pellet_mix':{},
              'pellet_phi':{}, 'pellet_r':{'length':1},
              'pellet_rate':{'particles':1,'time':-1},
              'pellet_rate_D2':{'particles':1,'time':-1},
              'pellet_ablrate':{'particles':1,'time':-1},
              'pellet_var':{'length':1}, 'pellet_var_tor':{'length':1},
              'pellet_velphi':{'velocity':1}, 'pellet_velr':{'velocity':1},
              'pellet_velz':{'velocity':1}, 'pellet_vx':{'velocity':1},
              'pellet_vy':{'velocity':1}, 'pellet_z':{'length':1},
              'r_p':{'length':1},
              }

    if custom is not None:
        expns.update(custom)
    elif trace in traces:
        expns.update(traces[trace])

    if units.lower()=='mks':
        if not isinstance(sim,fpy.sim_data):
            sim = fpy.sim_data(filename=filename)
        time   = unit_conv(trace_arr.time,   arr_dim='M3DC1', sim=sim, time=1)
        values = unit_conv(trace_arr.values, arr_dim='M3DC1', sim=sim, **expns)
    return fpy.sim_data.time_trace(values,time=time)



#-------------------------------------------
# Data type and array check ups and manipulation
#-------------------------------------------

# Find index of array elements with a specified value
def get_ind_at_val(arr, val, unique=True):
    #ToDo: Check type and length of val. In principle this can be a list
    ind = np.argwhere(arr==val)
    if unique:
        if len(ind[0])>1:
            raise Exception('Multiple occurences of '+str(val)+' found in '+str(arr)+'!')
        else:
            ind = ind[0,0]
    return ind


# Checks whether a value can be converted to a float
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

# Checks whether a value can be converted to an integer
def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False

def has_flux_coordinates(sim):
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename,time=time)
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    if (not isinstance(sim.fc,fpy.flux_coordinates)) or (fcoords!=None and (sim.fc.fcoords!=fcoords)) or (sim.fc.points!=points):
        if not fcoords:
            fcoords = ''
        sim = flux_coordinates(sim=sim, fcoords=fcoords, points=points, phit=phit,psin_range=psin_range)
    return sim



# Find value of array element that is nearest to a specified value
def find_nearest(arr, val):
    ind = np.abs(arr - val).argmin()
    return arr.flat[ind]


# Routines to check if a list is monotonic or strictly monotonic
def strict_inc(l):
    return all(x<y for x, y in zip(l, l[1:]))

def strict_dec(l):
    return all(x>y for x, y in zip(l, l[1:]))

def non_inc(l):
    return all(x>=y for x, y in zip(l, l[1:]))

def non_dec(l):
    return all(x<=y for x, y in zip(l, l[1:]))

def strict_monotonic(l):
    return strict_inc(l) or strict_dec(l)

def monotonic(l):
    return non_inc(l) or non_dec(l)



#-------------------------------------------
# Command line input / output routines
#-------------------------------------------

def prompt(message,options):
    """
    Read input from command line
    """
    input_str = ''
    if options == float:
        while not isinstance(input_str,float):
            input_str = input(message)
            if isfloat(input_str):
                input_str = float(input_str)
        #print('    '+str(input_str))
    elif options == int:
        while not isinstance(input_str,int):
            input_str = input(message)
            if isint(input_str):
                input_str = int(input_str)
    else:
        while not input_str in options:
            input_str = input(message)
        #print('    '+input_str)
    return input_str


# Routines to print certain types of messages (warning, error, note) to the command line
def printwarn(string):
    print(colored(string,'yellow'))
    return


def printerr(string):
    print(colored(string,'red'))
    return


def printnote(string):
    print(colored(string,'blue'))
    return


#-------------------------------------------
# File input / output. For reading hdf5 files please see read_h5.py.
#-------------------------------------------

def ReadTwoColFile(filename):
    """
    Reads a text file by column
    
    Arguments:

    **filename**
    Name of file to be read
    """
    with open(filename, 'r') as data:
        col1 = []
        col2 = []
        for line in data:
            r = line.split()
            col1.append(float(r[0]))
            col2.append(float(r[1]))
    col1 = np.asarray(col1)
    col2 = np.asarray(col2)
    return col1, col2



def ReadTwoColFile2(filename,header=0):
    """
    Reads a text file by column
    
    Arguments:

    **filename**
    Name of file to be read

    **header**
    Number of lines on top of file that are not part of the data columns
    """
    with open(filename, 'r') as data:
        col1 = []
        col2 = []
        if header > 0:
            colh = []
        i=0
        for line in data:
            r = line.split()
            if i < header:
                colh.append(r)
            else:
                
                if len(r)==2:
                    col1.append(float(r[0]))
                    col2.append(float(r[1]))
                else:
                    printwarn('WARNING: line omitted:'+line)
            i=i+1
    col1 = np.asarray(col1)
    col2 = np.asarray(col2)
    return colh, col1, col2



def get_filename(a):
    files = glob.glob(a)
    if len(files) < 1:
        raise Exception('No file found with name "'+a+'"!')
    else:
        if len(files) > 1:
            printwarn('WARNING: More than 1 file found. Using the newest one.')
            files.sort(key=os.path.getmtime,reverse=True)
        return files[0]


def get_lines(filename, linenumbers):
    """
    Reads and return only specific lines in a text file
    
    Arguments:

    **filename**
    Name of file to read

    **linenumbers**
    List of numbers of lines to read, e.g. [3,7,10], indexing starts at 0
    """
    lines = []
    with open(filename) as f:
        for i,line in enumerate(f):
            if i in linenumbers:
                lines.append(line)
                if len(lines)==len(linenumbers):
                    break
    return lines


def get_base_dirs(dirname):
    subfolders0 = [f.path for f in os.scandir(dirname) if (f.is_dir() and len(f.path)<6)]
    subfolders = []
    for dirname in list(subfolders0):
        subfolders.extend([f.path for f in os.scandir(dirname) if (f.is_dir() and 'base_' in f.path)])
    return subfolders

def get_run_dirs(dirname):
    subfolders0 = [f.path for f in os.scandir(dirname) if (f.is_dir() and len(f.path)<6)]
    subfolders = []
    for dirname in list(subfolders0):
        subfolders.extend([f.path for f in os.scandir(dirname) if (f.is_dir() and not 'base_' in f.path)])
    return subfolders
