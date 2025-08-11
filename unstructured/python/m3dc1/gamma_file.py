#!/usr/bin/env python3
import numpy as np

# Class that reads the text files containing information about growth rates from a toroidal mode number scan.
# Data includes some information about the physical model, equilibrium parameters and details about growth rate calculation

class Gamma_file:
    def __init__(self,filename):
        f = open(filename, 'r')
        data = f.readlines()
        f.close()
        datal1 = data[0].split()
        datal2 = data[1].split()
        
        for l,line in enumerate(data):
            #print(l,line)
            if 'gamma' in line:
                ind = l
                break
        
        width_height = False
        if ind == 3:
            width_height = True
        
        n_list = []
        gamma_list = []
        dgamma_list = []
        flat_list = []
        not_noisy_list = []
        gamma_manual_list = []
        gsconvgd = []
        finalerrgs = []
        pblist = []
        
        # Build lists from text file
        for i in range(ind+1,len(data)):
            n_list.append(int(data[i].split()[0]))
            gamma_list.append(float(data[i].split()[1])/2.0)
            dgamma_list.append(float(data[i].split()[2]))
            flat_list.append(int(data[i].split()[3]))
            not_noisy_list.append(int(data[i].split()[4]))
            gamma_manual_list.append(int(data[i].split()[5]))
            try:
                gsconvgd.append(int(data[i].split()[6]))
                finalerrgs.append(float(data[i].split()[7]))
            except:
                pass
            pblist.append(data[i].split()[8])
        
        self.vpnum = int(datal1[0])
        self.eta = float(datal1[1])
        self.bscale = float(datal1[2])
        self.rotation = int(datal1[3])
        self.fluidmodel = str(datal1[4])
        try:
            self.ipres = str(datal1[5])
        except:
            self.ipres = 0
            print('WARNING: ipres not found!')
        
        self.pped = float(datal2[0])
        self.jped = float(datal2[1])
        if len(datal2)>=3:
            self.jeliteped = float(datal2[2])
        else:
            print('WARNING: jelite not found!')
        if len(datal2)>=4:
            self.omegsti_max = float(datal2[3]) #Diamagnetic frequency NOT mulitplied by n, i.e. d pi / d psi / (ei ni)
        else:
            self.omegsti_max = -1
            print('WARNING: diamagnetic frequency not found!')
        if ind > 2:
            datal3 = data[2].split()
            self.ped_height = datal3[0]
            self.ped_width = datal3[1]
        self.n_list = np.asarray(n_list,dtype=int)
        self.gamma_list = np.asarray(gamma_list)
        if self.omegsti_max > 0:
            self.gamma_dia_list = self.gamma_list / (self.n_list * self.omegsti_max/4)
        self.dgamma_list = np.asarray(dgamma_list)
        self.flat_list = np.asarray(flat_list,dtype=int)
        self.not_noisy_list = np.asarray(not_noisy_list,dtype=int)
        self.gamma_manual_list = np.asarray(gamma_manual_list,dtype=int)
        self.gsconvgd = np.asarray(gsconvgd,dtype=int)
        self.finalerrgs = np.asarray(finalerrgs)
        self.pblist = np.asarray(pblist,dtype=int)
