import eccentricdark as ed
import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate 


def generate_invcdf(pdf, xmin, xmax, xscale): 

    #calculate cdf
    if xscale=='log':
        xvals = np.logspace(np.log10(xmin), np.log10(xmax), ed.invcdf_density_factor)
        xmidpoints = (xvals[1:] + xvals[:-1]) / 2.
        dxvals = np.diff(xvals)

        normalization = np.sum(pdf(xmidpoints)*dxvals)
        normed_pdf = lambda x: pdf(x)/normalization   

        cdf_vals = np.zeros(len(xvals)+1)
        cdf_vals[1:-1] = np.cumsum(normed_pdf(xmidpoints)*dxvals)
        cdf_vals[-1] = 1.


        fullxvals = np.append(np.append(xvals[0], xmidpoints), xvals[-1])

        cdf_interp = scipy.interpolate.interp1d(fullxvals, cdf_vals)

        log10invcdf_interp = scipy.interpolate.interp1d(cdf_vals, fullxvals)
        
        invcdf_interp = lambda x: log10invcdf_interp(x)

    elif xscale=='linear': 
        normalization = scipy.integrate.quad(pdf, xmin, xmax)[0]
        normed_pdf = lambda x: pdf(x)/normalization
        #size_xbasis = np.int(np.floor((xmax-xmin)*ed.invcdf_density_factor))
        xeval = np.linspace(xmin, xmax, ed.invcdf_density_factor) 

        cdf_func = np.vectorize(lambda x: scipy.integrate.quad(normed_pdf, xmin, x)[0])
        cdf = cdf_func(xeval)
        cdf_interp = scipy.interpolate.interp1d(xeval, cdf) 

        invcdf_interp = scipy.interpolate.interp1d(cdf, xeval) 

    elif xscale == 'test':
        normalization = scipy.integrate.quad(pdf, xmin, xmax)[0]
        normed_pdf = lambda x: pdf(x)/normalization

        xvals = np.logspace(np.log10(xmin), np.log10(xmax), ed.invcdf_density_factor)
        cdf = lambda x: scipy.integrate.quad(normed_pdf, xmin, x)[0]
        cdf_vals = np.vectorize(cdf)(xvals)
        
        inv_cdf = scipy.interpolate.interp1d(cdf_vals, xvals)
        invcdf_interp = lambda x: inv_cdf(x)
        
        plt.plot(np.log10(xvals), normed_pdf(xvals), label='True PDF')


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





