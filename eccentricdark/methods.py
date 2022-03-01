import eccentricdark as ed
import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate 


def generate_invcdf(pdf, xmin, xmax): 
    #normalize pdf
    normalization = scipy.integrate.quad(pdf, xmin, xmax)[0] 
    normed_pdf = lambda x: pdf(x)/normalization

    #calculate cdf
    if ((np.log10(xmax)-np.log10(xmin))>3): 
        xscale='log'
        size_xbasis = np.int(np.floor((np.log10(xmax)-np.log10(xmin))*ed.invcdf_density_factor))
        xeval = np.logspace(np.log10(xmin), np.log10(xmax), size_xbasis)
    else: 
        xscale='linear'
        size_xbasis = np.int(np.floor((xmax-xmin)*ed.invcdf_density_factor))
        xeval = np.linspace(xmin, xmax, size_xbasis) 

    cdf_func = np.vectorize(lambda y: scipy.integrate.quad(normed_pdf, xmin, y)[0])
    cdf = cdf_func(xeval)
    cdf_interp = scipy.interpolate.interp1d(xeval, cdf) 
        #CAREFUL - might want to interpolate log of values for log spacing

    #invert cdf 
    invcdf_interp = scipy.interpolate.interp1d(cdf, xeval) 
        #CAREFUL - might want to interpolate log of values for log spacing

    return invcdf_interp

def multigauss(y, u, a, b, c, k): 
    return (
        (k/y)
        * np.exp(
            a * np.power(-u + np.log(y), 2.)
            + b * np.power(-u + np.log(y), 3.)
            + c * np.power(-u + np.log(y), 4.)
        )
    )





