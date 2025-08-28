#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 17 2025

@author: Andreas Kleiner
"""
import numpy as np
import os.path

def find_str_in_array(s,array,return_linenum=False):
    """
    Find a string in an array that consists of strings
    and return the first element that contains this string.
    Optionally, also return array index (line number).

    Arguments:

    **s**
    String to be searched for.

    **array**
    Array consisting of strings, e.g. a text file with one string
    representing each line.

    **return_linenum**
    True/False. If True, return array index where string was found in
    addition to the string.
    """
    string=''
    for n,line in enumerate(array):
        if s in line and line[0] != '!':
            string = line
            linenum = n
            break
    if return_linenum:
        return string, linenum
    else:
        return string


def convert_meshgen_input(meshgen_filename='input',mfmgen_filename='input_mfmgen'):
    """
    Converts input file for m3dc1_meshgen to format used by m3dc1_mfmgen. The converted
    input file is being written to a new file.
    Note: Works for modelType 0 and 3. Other model types can be implemented, if required.

    Arguments:

    **meshgen_filename**
    Path and filename to meshgen input file.

    **mfmgen_filename**
    Filename for m3dc1_mfmgen input file.
    """
    
    #Read m3dc1_meshgen input file
    f = open(meshgen_filename, 'r')
    infile = f.readlines()
    f.close()
    
    search_for = ['modelType','inFile','outFile','meshSize','useVacuumParams','adjustVacuumParams','vacuumParams','numVacuumPts','meshGradationRate','resistive-width','plasma-offsetX','plasma-offsetY','vacuum-width','vacuum-height']
    
    input_values = {}
    
    for s in search_for:
        value = find_str_in_array(s,infile).replace('\n','')
        input_values[s] = value
        #print(value)
    
    print('The following values were read from the input file:')
    print(input_values)
    print('\n\n')
    
    modelType = input_values['modelType'].split(' ')[1]
    if modelType not in ['3','0']:
        print('ERROR: modelType not recognized. Please check input file!')
        return
    #print(modelType)
    output_lines = []
    
    # Read quantities that are used by both modelType 0 and 3
    outFile = input_values['outFile'].split(' ')[1]
    meshSize = input_values['meshSize'].split(' ')[1:]
    vacuumParams = input_values['vacuumParams'].split(' ')[1:]
    numVacuumPts = input_values['numVacuumPts'].split(' ')[1]
    
    if modelType == '3':
        inFile = input_values['inFile'].split(' ')[1]
        useVacuumParams = input_values['useVacuumParams'].split(' ')[1]
        adjustVacuumParams = input_values['adjustVacuumParams'].split(' ')[1]
        meshGradationRate = input_values['meshGradationRate'].split(' ')[1]
        resistivewidth = input_values['resistive-width'].split(' ')[1]
        plasmaoffsetX = input_values['plasma-offsetX'].split(' ')[1]
        plasmaoffsetY = input_values['plasma-offsetY'].split(' ')[1]
        vacuumwidth = input_values['vacuum-width'].split(' ')[1]
        vacuumheight = input_values['vacuum-height'].split(' ')[1]
        
        #print(meshSize)
        
        #Write m3dc1_mfmgen file
        output_lines.append("numBdry 1\n")
        output_lines.append("bdryFile " + inFile + " 1 "+str(meshSize[0]) +"\n\n")
        output_lines.append("thickWall 1 1 2 "+str(resistivewidth)+"\n\n")
        output_lines.append("useVacuum 2 3 0.1\n")#ToDo: check values
        output_lines.append("vacuumParams " + ' '.join(vacuumParams) + "\n")
        output_lines.append("numVacuumPts " + str(numVacuumPts) + "\n\n")
        output_lines.append("numFace 3\n")
        output_lines.append("faceBdry  1 1 "+str(meshSize[0])+"\n")
        output_lines.append("faceBdry  2 1 2 "+str(meshSize[1])+"\n")
        output_lines.append("faceBdry  2 2 3 "+str(meshSize[2])+"\n\n")
        
        output_lines.append("meshGradationRate " + str(meshGradationRate) + "\n\n")
        output_lines.append("outFile " + outFile + "\n")
        
        print("Make sure the following lines have been added to C1input:\n\n"+ \
            "\tboundary_type(1) = 1  ! First wall inner boundary\n"+ \
            "\tboundary_type(2) = 0  ! First wall outer boundary\n"+ \
            "\tboundary_type(3) = 2  ! Domain boundary\n"+ \
            "\tzone_type(1) = 1      ! Plasma\n"+ \
            "\tzone_type(2) = 2      ! First wall\n"+ \
            "\tzone_type(3) = 3      ! Vacuum")
    elif modelType == '0':
        #Write m3dc1_mfmgen file
        output_lines.append("numBdry 0\n\n")
        output_lines.append("useVacuum 2 1 " + str(meshSize[0]) + "\n")
        output_lines.append("vacuumParams " + ' '.join(vacuumParams) + "\n")
        output_lines.append("numVacuumPts " + str(numVacuumPts) + "\n\n")
        output_lines.append("numFace 1\n")
        output_lines.append("faceBdry  1 1 "+str(meshSize[0])+"\n\n")
        output_lines.append("outFile " + outFile + "\n")
        
        print("Make sure the following lines have been added to C1input:\n\n"+ \
            "\tboundary_type(1) = 2  ! Domain boundary\n"+ \
            "\tzone_type(1) = 1      ! Plasma")
    
    #print(output_lines)
    
    if not os.path.isfile(mfmgen_filename):
        f = open(mfmgen_filename, 'w')
        for line in output_lines:
            f.write(line)
        f.close()
    else:
        print('ERROR: File ' + mfmgen_filename + ' already exists. Please choose another filename!')
    return

