#!/usr/bin/env python3
import numpy as np

# This class contains information about the growth rates calculated from multiple linear simulations with varying toroidal mode number.
class Gamma_data:
    def __init__(self,n_list, gamma_list, dgamma_list, flat_list, not_noisy_list, gamma_manual_list, gsconvgd, finalerrgs, pblist, ped_loc):
        
        self.n_list = np.asarray(n_list)
        self.gamma_list = np.asarray(gamma_list)
        self.dgamma_list = np.asarray(dgamma_list)
        self.flat_list = np.asarray(flat_list)
        self.not_noisy_list = np.asarray(not_noisy_list)
        self.gamma_manual_list = np.asarray(gamma_manual_list)
        self.gsconvgd = np.asarray(gsconvgd)
        self.finalerrgs = np.asarray(finalerrgs)
        self.pblist = np.asarray(pblist)
        self.ped_loc = np.asarray(ped_loc)
