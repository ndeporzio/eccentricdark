# The functions here define the binary eccentricity distribution functions...

import scipy
import scipy.stats
import numpy as np
import eccentricdark as ed

def estar_sampler(
    estar_distribution,
    args=None 
): 
    #should return the invcdf function for specified pdf 

    if estar_distribution=="fixed": 
        return np.vectorize(lambda x: args) 

    if estar_distribution=="loggaussian":
        # args[0] = center of gaussian in log space
        # args[1] = width of gaussian in log space
        return np.vectorize(
            lambda x: np.power(10., scipy.stats.norm(loc=args[0], scale=args[1]).ppf(x)))

    if estar_distribution=="isolated": 
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar, 
            u = -12.308577637842925, 
            a = -0.3000192097507729,
            b = -0.02842041062114348, 
            c = -0.0024862991498544496, 
            k = 0.23232375892104137
        ) 

        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-2))

    if estar_distribution=="ejected": 
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar, 
            u = -8.903606734106413, 
            a = -0.21505957280762572,
            b = -0.05215442618559087, 
            c = -0.005523661062693592, 
            k = 0.2504431992363881
        ) 

        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-2))

    if estar_distribution=="incluster": 
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar, 
            u = -7.1814548324724425, 
            a = -0.41160569016079956,
            b = -0.13633206068723644, 
            c = -0.01848512571848316, 
            k = 0.34149918878761343
        ) 


        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-1))

    if estar_distribution=="galcenter": 
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar, 
            u = -8.670877869456227, 
            a = -1.4033086813787914,
            b = -1.0852555343749364, 
            c = -0.3468739524658422, 
            k = 0.06114451273864072
        ) 


        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-1))
