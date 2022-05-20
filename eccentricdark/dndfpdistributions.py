import scipy
import scipy.stats
import numpy as np
import eccentricdark as ed

def dndfp_sampler(
    mc, 
    e_of_fp, 
    fpmin, 
    fpmax, 
    xscale='log'
):

    pdf = (lambda fp: np.power(fp, -11./3.) * ed.F_script(e_of_fp(fp)))

    return ed.generate_invcdf(pdf, fpmin, fpmax, xscale)(np.random.random())

    

