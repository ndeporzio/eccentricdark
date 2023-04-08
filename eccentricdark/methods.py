import eccentricdark as ed
import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate 


def generate_invcdf(pdf, xmin, xmax, xscale): 

    if xscale=='log':
         xvals = np.linspace(np.log10(xmin), np.log10(xmax), ed.invcdf_density_factor)
         xmidpoints = (xvals[1:] + xvals[:-1]) / 2.
         dxvals = np.diff(xvals)

         normalization = np.sum(pdf(np.power(10., xmidpoints))*dxvals)
         normed_pdf = lambda x: pdf(x)/normalization   

         cdf_vals = np.zeros(len(xvals)+1)
         cdf_vals[1:-1] = np.cumsum(normed_pdf(np.power(10., xmidpoints))*dxvals)
         cdf_vals[-1] = 1.

         fullxvals = np.append(np.append(xvals[0], xmidpoints), xvals[-1])

         cdf_interp = scipy.interpolate.interp1d(fullxvals, cdf_vals)
         log10invcdf_interp = scipy.interpolate.interp1d(cdf_vals, fullxvals)
         invcdf_interp = lambda x: np.power(10., log10invcdf_interp(x))

        ########################

    elif xscale=='linear': 
        normalization = scipy.integrate.quad(pdf, xmin, xmax)[0]
        normed_pdf = lambda x: pdf(x)/normalization

        xeval = np.linspace(xmin, xmax, 100)

        cdf_func = np.vectorize(lambda x: scipy.integrate.quad(normed_pdf, xmin, x)[0])
        cdf = cdf_func(xeval)
        
        deleteidxs = []
        for idx, val in enumerate(cdf[:-1]): 
            if ((val==cdf[idx+1]) and (val==0.)): 
                deleteidxs.append(idx)
            elif ((val==cdf[idx+1]) and (val==1.)):
                deleteidxs.append(idx+1)
                
        xeval = np.delete(xeval, deleteidxs)
        cdf = np.delete(cdf, deleteidxs)

        cdf_interp = scipy.interpolate.interp1d(xeval, cdf)
        invcdf_interp = scipy.interpolate.interp1d(cdf, xeval)

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





