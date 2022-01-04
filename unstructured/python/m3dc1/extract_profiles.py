import numpy as np
import os


def extract_profiles(filename):
    os.system('extract_profiles.sh '+filename)
    
    if os.path.isfile('profile_p'):
        convert_p()
        print('profile_p converted to Pascal')
    
    return


def convert_p():
    
    f = open('profile_p', 'r')
    data = f.readlines()
    f.close()
    nlines = len(data)
    
    profile = np.empty((nlines,3), dtype=float)
    
    for i,line in enumerate(data):
        profile[i,:] = list(map(float, line.split()))
    profile[:,1] = profile[:,1]*1000
    profile[:,2] = profile[:,2]*1000
    print(profile)
    
    with open('profile_p', 'w') as pf:
        for ln in profile:
            pf.write(" {0:8.6f}   {1:>9.3f}   {2:>9.3f}\n".format(*ln))
    
    return
