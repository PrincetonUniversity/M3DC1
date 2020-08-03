#!/usr/bin/env python3
import numpy as np


class Gamma_file:
    def __init__(self,filename):
        f = open(filename, 'r')
        data = f.readlines()
        datal1 = data[0].split()
        datal2 = data[1].split()
        
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
        for i in range(3,len(data)):
            n_list.append(int(data[i].split()[0]))
            gamma_list.append(float(data[i].split()[1]))
            dgamma_list.append(float(data[i].split()[2]))
            flat_list.append(int(data[i].split()[3]))
            not_noisy_list.append(int(data[i].split()[4]))
            gamma_manual_list.append(int(data[i].split()[5]))
            try:
                gsconvgd.append(int(data[i].split()[6]))
                finalerrgs.append(float(data[i].split()[7]))
            except:
                pass
            pblist.append(int(data[i].split()[8]))
        
        self.vpnum = int(datal1[0])
        self.eta = float(datal1[1])
        self.bscale = float(datal1[2])
        self.rotation = int(datal1[3])
        self.fluidmodel = str(datal1[4])
        self.jped = float(datal2[0])
        self.pped = float(datal2[1])
        self.n_list = np.asarray(n_list)
        self.gamma_list = np.asarray(gamma_list)
        self.dgamma_list = np.asarray(dgamma_list)
        self.flat_list = np.asarray(flat_list)
        self.not_noisy_list = np.asarray(not_noisy_list)
        self.gamma_manual_list = np.asarray(gamma_manual_list)
        self.gsconvgd = np.asarray(gsconvgd)
        self.finalerrgs = np.asarray(finalerrgs)
        self.pblist = np.asarray(pblist)
