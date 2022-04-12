# The functions here define the binary eccentricity distribution functions...

import scipy
import numpy as np
import eccentricdark as ed

def chi_distribution_sampler(
    form,
    args=None
):

    if (form=="fixed"):
        chi = args[0] 
        return chi 
    elif (form=="uniform"):
        chi_min = 0.
        chi_max = args[0] 
        return chi_max*np.power(np.random.random(), 1./3.)

    

