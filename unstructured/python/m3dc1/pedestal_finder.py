import numpy as np
import os
import matplotlib.pyplot as plt
import m3dc1.fpylib as fpyl
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)


def Hmode_profiles(edge=0.08, ped=0.4, core=2.5, rgrid=201, expin=1.5, expout=1.5, widthp=0.04, fac=1.0, xphalf=None):
    """
     This function generates H-mode  density and temperature profiles evenly
     spaced in your favorite radial coordinate

    :param edge: (float) separatrix height

    :param ped: (float) pedestal height

    :param core: (float) on-axis profile height

    :param rgrid: (int) number of radial grid pointsx

    :param expin: (float) inner core exponent for H-mode pedestal profile

    :param expout (float) outer core exponent for H-mode pedestal profile

    :param width: (float) width of pedestal

    :param xphalf: (float) position of tanh
    """

    w_E1 = 0.5 * widthp  # width as defined in eped
    if xphalf is None:
        xphalf = 1.0 - fac*w_E1
    xphalf0 = 1.0 - w_E1

    xped = xphalf - w_E1

    pconst = 1.0 - np.tanh((1.0 - xphalf0) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + np.tanh(1.0) - pconst)

    coretanh = 0.5 * a_t * (1.0 - np.tanh(-xphalf / w_E1) - pconst) + edge

    xpsi = np.linspace(0, 1, rgrid)
    ones = np.ones(rgrid)

    val = 0.5 * a_t * (1.0 - np.tanh((xpsi - xphalf) / w_E1) - pconst) + edge * ones

    #val = 1.0 - np.tanh((xpsi - xphalf) / w_E1)
    xtoped = xpsi / xped
    for i in range(0, rgrid):
        if xtoped[i] ** expin < 1.0:
            val[i] = val[i] + (core - coretanh) * (1.0 - xtoped[i] ** expin) ** expout

    return val


def Hmode_profiles_fit(xpsi,edge=0.08, ped=0.4, core=2.5, expin=1.5, expout=1.5, widthp=0.04):
    """
     This function generates H-mode  density and temperature profiles evenly
     spaced in your favorite radial coordinate

    :param edge: (float) separatrix height

    :param ped: (float) pedestal height

    :param core: (float) on-axis profile height

    :param rgrid: (int) number of radial grid pointsx

    :param expin: (float) inner core exponent for H-mode pedestal profile

    :param expout (float) outer core exponent for H-mode pedestal profile

    :param width: (float) width of pedestal

    :param xphalf: (float) position of tanh
    """
    fac=1.0
    rgrid=len(xpsi)
    xphalf=None
    w_E1 = 0.5 * widthp  # width as defined in eped
    if xphalf is None:
        xphalf = 1.0 - fac*w_E1
    xphalf0 = 1.0 - w_E1

    xped = xphalf - w_E1

    pconst = 1.0 - np.tanh((1.0 - xphalf0) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + np.tanh(1.0) - pconst)

    coretanh = 0.5 * a_t * (1.0 - np.tanh(-xphalf / w_E1) - pconst) + edge

    #xpsi = np.linspace(0, 1, rgrid)
    ones = np.ones(rgrid)

    val = 0.5 * a_t * (1.0 - np.tanh((xpsi - xphalf) / w_E1) - pconst) + edge * ones

    #val = 1.0 - np.tanh((xpsi - xphalf) / w_E1)
    xtoped = xpsi / xped
    for i in range(0, rgrid):
        if xtoped[i] ** expin < 1.0:
            val[i] = val[i] + (core - coretanh) * (1.0 - xtoped[i] ** expin) ** expout

    return val

def pedestal_finder(profile, psi_norm=None, eped_definition=False, return_fit=False, doPlot=False, ngrid=201):
    from scipy import optimize, interpolate

    # param : plasma profile,  a 1D profile on a rho or psi grid (i.e. Te, P, ni)
    # param : return_fit, returns profile_fit_grid and profile_fit
    # param : doPlot, compares the fit to the data

    # if return_fit=False: return pedestal top, pedestal width
    # if return_fit=True:  return pedestal top, pedestal width, profile_fit_grid, profile_fit
    if isinstance(profile, int):
        printe('ERROR! ne input should be a 1D profile')

    def cost_function(c):
        if any(c < 0):
            return 1e10

        if eped_definition:
            nval = Hmode_profiles(
                rgrid=len(profile), xphalf=None, widthp=c[1], core=profile[0], ped=c[0], edge=profile[-1], expin=2.0, expout=2.0
            )
        else:
            nval = Hmode_profiles(
                rgrid=len(profile), xphalf=c[2], widthp=c[1], core=profile[0], ped=c[0], edge=profile[-1], expin=2.0, expout=2.0
            )

        cost = np.sqrt(sum(((profile - nval) ** 2 / profile[0] ** 2) * weight_func))
        return cost

    if psi_norm is None:
        psi_norm = np.linspace(0, 1, len(profile))

    profile = np.interp(fp=profile, xp=psi_norm, x=np.linspace(0, 1, ngrid))
    psi_norm = np.linspace(0, 1, ngrid)

    inversion_point = np.argmin(np.gradient(profile[int(ngrid / 1.6) :])) + int(ngrid / 1.6)
    inversion_point_margin = inversion_point - int(0.1 * ngrid)

    weight_func = np.zeros(ngrid)
    weight_func[inversion_point_margin:] += 1.0

    width0 = 1 - psi_norm[inversion_point]

    # eped definition sets xphalf to 1 - width
    # non eped definition optimizes over xphalf
    if eped_definition:
        c = [interpolate.interp1d(psi_norm, profile)([1 - 2 * width0])[0], width0]
        c = list(map(float, optimize.minimize(cost_function, c, method='Nelder-Mead', jac=False).x))
        nval_fit = Hmode_profiles(rgrid=ngrid, widthp=c[1], core=profile[0], ped=c[0], xphalf=None, edge=profile[-1], expin=2.0, expout=2.0)
    else:
        xphalf0 = 1 - width0 * 0.5
        c = [interpolate.interp1d(psi_norm, profile)([1 - 2 * width0])[0], width0, xphalf0]
        c = list(map(float, optimize.minimize(cost_function, c, method='Nelder-Mead', jac=False).x))
        nval_fit = Hmode_profiles(rgrid=ngrid, widthp=c[1], core=profile[0], ped=c[0], xphalf=c[2], edge=profile[-1], expin=2.0, expout=2.0)

    if doPlot:
        plt.figure()
        plt.plot(psi_norm, profile, label='raw')
        plt.plot(psi_norm, nval_fit, label='fit')
        plt.legend()

    if return_fit:
        return c[0], c[1], psi_norm, nval_fit  # pedestal top, pedestal width, profile_fit_grid, profile_fit
    return c[0], c[1]  # pedestal top, pedestal width



def pedestal_finder_fit(profile, psi_norm=None,psin_cutoff=0.7,doPlot=True):
    psin_min_ind = fpyl.get_ind_near_val(psi_norm,psin_cutoff)
    sep_height = profile[-1]
    ped_height = profile[psin_min_ind]
    core_height = profile[0]
    
    print(sep_height,ped_height,core_height)
    print('fitting with psin_cutoff='+ str(psin_cutoff))
    fit = curve_fit(Hmode_profiles_fit,psi_norm[psin_min_ind:],profile[psin_min_ind:],[sep_height,ped_height,core_height,1.1,1.9,0.2])
    temp = Hmode_profiles(fit[0][0],fit[0][1],fit[0][2],len(psi_norm),fit[0][3],fit[0][4],fit[0][5],1.0,None)
    #print(fit)
    
    if doPlot:
        plt.figure()
        plt.plot(psi_norm,profile)
        plt.plot(np.linspace(0, 1, len(psi_norm)),temp)
        plt.title('/'.join(os.getcwd().split('/')[-3:]))
    return fit[0][1],fit[0][5]



def get_ped_structure(profile, psi_norm=None,fit=True,psin_cutoff=0.7,ngrid=200,doPlot=False):
    if fit:
        sign = 1
        psin_cutoff_init = psin_cutoff
        while True:
            if psin_cutoff<0 or psin_cutoff>1:
                fpyl.printerr('ERROR: psin out of range. Fit did not converge!')
                break
            try:
                ped_top,ped_wid = pedestal_finder_fit(profile,psi_norm=psi_norm,psin_cutoff=psin_cutoff,doPlot=doPlot)
                break
            except (RuntimeError,OptimizeWarning):
                fpyl.printerr('Error: Fit did not converge!')
            #First increase psin_cutoff until it reaches 1, then decrease from initial value
            if psin_cutoff >= 0.87:
                print('Decreasing cutoff for pedestal fit...')
                sign = -1
                psin_cutoff = psin_cutoff_init
            
            psin_cutoff = psin_cutoff + sign*0.05
        
    else:
        ped_top,ped_wid = pedestal_finder(profile,psi_norm=psi_norm,ngrid=ngrid)
    return ped_top,ped_wid
