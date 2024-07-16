#!/usr/bin/env python3
import numpy as np


def compensate_renorm(x):
    """
    Compensate for energy renormalization in linear simulations.

    Arguments:

    **x**
    The array to be renormalized.
    """
    y = np.asarray(x)
    y = y.astype(np.float64)
    
    for k in range(len(x)-2):
        if not (x[k]*x[k+1] <= 0.):
            if (abs(x[k+1]/x[k]) < 1e-9):
                print('Renormalization found at '+str(k))
                # Rescale
                if(k < len(x)-2):
                    dx = x[k+2] - x[k+1]
                else:
                    dx = 0.
                # this is what x[k] would be assuming exponential growth
                #print(dx,x[k+1],-dx/(x[k+1]+dx/2.))
                f = x[k+1]*np.exp(-dx/(x[k+1]+dx/2.))
                #print(f)
                for j in range(k+1):
                    #print(y[j],f, x[k])
                    y[j] = y[j] * f 
                    y[j] = y[j] / x[k]
                    #print(y[j])
    return y

