import fpy
import matplotlib.pyplot as plt
from unit_conv import unit_conv


def plot_diagnostics(file_name='C1.h5',units='mks'):
    """
    This plots three diagnostics.
    1 - The time needed per timestep
    2 - The number of iterations needed per timestep
    3 - The temporal values of each time-slice
    """
    
    sim     =  fpy.sim_data(file_name)
    t_slices=  sim.get_diagnostic('slice times')
    timings =  sim.get_diagnostic('timings')
    kspits  =  sim.get_diagnostic('iterations')

    
    
    f = plt.figure(1)
    plt.subplot(121)
    sum = []
    for string in list(timings.diagnostic.keys()):
        if any(timings.diagnostic[string].value) != 0:
            string_label = string.replace('_', '\_')
            plt.plot(timings.x_axis, timings.diagnostic[string].value, label='{}'.format(string_label), alpha=0.5)
        if len(sum) == 0:
            sum = timings.diagnostic[string].value
        sum += timings.diagnostic[string].value
    plt.plot(timings.x_axis,sum,label='Total time')
    plt.xlabel('timestep')
    plt.ylabel('Time taken [s]')
    plt.legend(loc='best')
    
    
    plt.subplot(122)
    plt.plot(kspits.x_axis, kspits.diagnostic,)
    plt.ylabel('Number of iterations')
    plt.xlabel('timestep')

    f.show()


    g = plt.figure(2)
    x_axis = t_slices.x_axis
    plt.xlabel('Time')
    if units=='mks':
        x_axis = unit_conv(x_axis,time=1)
        plt.xlabel('Time [s]')
    plt.scatter(x_axis,t_slices.diagnostic)
    plt.ylabel('Time slice label')
    g.show()



    
