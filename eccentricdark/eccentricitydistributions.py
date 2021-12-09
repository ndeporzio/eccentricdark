# The functions here define the binary eccentricity distribution functions...

import scipy
import numpy as np
import eccentricdark as ed

def estar_sampler(
    estar_distribution,
    args=None
): 

    def multigauss(y, u, a, b, c, k): 
        return (
            (k/y)
            * np.exp(
                a * np.power(-u + np.log(y), 2.)
                + b * np.power(-u + np.log(y), 3.)
                + c * np.power(-u + np.log(y), 4.)
            )
        )

    if estar_distribution=="fixed": 
        estar_sampled = args
        return estar_sampled

    if estar_distribution=="gaussian":
        # args[0] = center of gaussian in log space
        # args[1] = width of gaussian in log space
        return np.power(10., np.random.normal(args[0], args[1]))

    if estar_distribution=="isolated": 
        pdf_unnormed = lambda estar : multigauss(
            y = estar, 
            u = -12.3086, 
            a = -0.300019,
            b = -0.0284204, 
            c = -0.0024863, 
            k = 0.232324
        ) 

        pdf_normed = lambda estar: (
            pdf_unnormed(estar)
            / scipy.integrate.quad(pdf_unnormed, 1.0e-12, 1.0e-1)
        ) 
        
        random_sample = np.random.random(1)
        
        counter = 0.
        estar = np.power(10., -12)
        while (counter < random_sample): 
            counter += np.power(10., 0.1)*pdf_normed(estar) 

        return estar 

