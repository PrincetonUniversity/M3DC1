#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: akleiner
"""
import numpy as np
import math
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
    if equal==True:
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
    
    if (type(x) is np.ndarray or type(x) is list):
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



#-------------------------------------------
# Unit conversion
#-------------------------------------------


# Returns field label depending on chosen system of units
def get_fieldlabel(units,field):
    if units.lower()=='m3dc1':
        if field == 'j':
            label = 'current density'
        elif field == 'ni':
            label = 'ion density'
        elif field == 'ne':
            label = 'electron density'
        elif field == 'v':
            label = 'velocity'
        elif field == 'B':
            label = 'magnetic field strength'
        elif field == 'p':
            label = 'pressure'
        elif field == 'pi':
            label = 'ion pressure'
        elif field == 'pe':
            label = 'electron pressure'
        elif field == 'ti':
            label = 'ion temperature'
        elif field == 'te':
            label = 'electron temperature'
        elif field == 'A':
            label = 'vector potential'
        elif field == 'gradA':
            label = 'grad vector potential'
        elif field == 'E':
            label = 'electric field '
        else:
            label = field
        unit_label = 'M3DC1 units'
    
    if units.lower()=='mks':
        if field == 'j':
            label = 'current density'
            unit_label = '$A/m^2$'
        elif field == 'ni':
            label = 'ion density'
            unit_label = 'particles/$m^3$'
        elif field == 'ne':
            label = 'electron density'
            unit_label = 'particles/$m^3$'
        elif field == 'v':
            label = 'velocity'
            unit_label = '$m/s$'
        elif field == 'B':
            label = 'magnetic field strength'
            unit_label = '$T$'
        elif field == 'p':
            label = 'pressure'
            unit_label = '$Pa$'
        elif field == 'pi':
            label = 'ion pressure'
            unit_label = '$Pa$'
        elif field == 'pe':
            label = 'electron pressure'
            unit_label = '$Pa$'
        elif field == 'ti':
            label = 'ion temperature'
            unit_label = '$eV$'
        elif field == 'te':
            label = 'electron temperature'
            unit_label = '$eV$'
        elif field == 'A':
            label = 'vector potential'
            unit_label = '$Tesla \cdot m$'
        elif field == 'gradA':
            label = 'grad vector potential'
            unit_label = '$(Tesla \cdot m)$ / (m or rad)'
        elif field == 'E':
            label = 'electric field'
            unit_label = '$V/m$'
        else:
            label = field
            unit_label = 'MKS units'
        #label = field + ' (' + label + ')'
    return label, unit_label


def get_conv_field(units,field,field1_ave):
    """
    Returns converted field depending on chosen system of units
    """
    if units.lower()=='m3dc1':
        if field == 'j':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',current_density=1)
        elif field == 'ni':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',particles=1,length=-3)
        elif field == 'ne':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',particles=1,length=-3)
        elif field == 'v':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',velocity=1)
        elif field == 'B':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',magnetic_field=1)
        elif field == 'p':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
        elif field == 'pi':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
        elif field == 'pe':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
        elif field == 'ti':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',temperature=1)
        elif field == 'te':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',temperature=1)
        elif field == 'A':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',magnetic_field=1,length=1)
        elif field == 'grad A':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',magnetic_field=1)#ToDo: Normalization for derivative wrt phi
        elif field == 'E':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',electric_field=1)
    return field1_ave


def get_conv_trace(units,trace,trace_arr):
    """
    Returns converted time trace depending on chosen system of units
    """
    #ToDo: Add missing traces
    if units.lower()=='mks':
        if trace == 'time':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',time=1)
        elif trace == 'ke':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',energy=1)
        elif trace == 'me':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',energy=1)
        elif trace == 'p':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',energy=1)
        elif trace == 'ip':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',current=1)
        elif trace == 'it':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',current=1)
        elif trace == 'iw':
            trace_arr = unit_conv(trace_arr,arr_dim='M3DC1',current=1)
        else:
            printwarn('WARNING: Unit conversion not yet implemented for ' + trace + '. Feel free to add it.')
    return trace_arr



#-------------------------------------------
# Data type and array check ups and manipulation
#-------------------------------------------

# Find index of array elements with a specified value
def get_ind_at_val(arr, val, unique=True):
    #ToDo: Check type and length of val. In principle this can be a list
    ind = np.argwhere(arr==val)
    if unique==True:
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


# Find index of array element that is nearest to a specified value
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
        while type(input_str) is not float:
            input_str = input(message)
            if isfloat(input_str)==True:
                input_str = float(input_str)
        #print('    '+str(input_str))
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

def ReadTwoColFile(file_name):
    """
    Reads a text file by column
    
    Arguments:

    **file_name**
    Name of file to be read
    """
    with open(file_name, 'r') as data:
        col1 = []
        col2 = []
        for line in data:
            r = line.split()
            col1.append(float(r[0]))
            col2.append(float(r[1]))
    col1 = np.asarray(col1)
    col2 = np.asarray(col2)
    return col1, col2



def ReadTwoColFile2(file_name,header=0):
    """
    Reads a text file by column
    
    Arguments:

    **file_name**
    Name of file to be read

    **header**
    Number of lines on top of file that are not part of the data columns
    """
    with open(file_name, 'r') as data:
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
