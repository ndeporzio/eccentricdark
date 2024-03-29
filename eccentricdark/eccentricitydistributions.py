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
        # args = value distribution is fixed to 
        return np.vectorize(lambda x: args) 

    if estar_distribution=="loggaussian":
        # args[0] = center of gaussian in log space
        # args[1] = width of gaussian in log space
        return np.vectorize(
            lambda x: np.power(10., scipy.stats.norm(loc=args[0], scale=args[1]).ppf(x)))

    if estar_distribution=="isolated": 
        fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss, 
            np.power(10., ed.fieldData[:,0]),
            ed.fieldData[:,1]/max(ed.fieldData[:,1]),
            [-12.4, -0.43, -0.10, -1.4e-2, 1.6e-6]
        )
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        )
        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-2, 'log'))

    if estar_distribution=="ejected": 
        fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss, 
            np.power(10., ed.ejectedData[:,0]), 
            ed.ejectedData[:,1]/max(ed.ejectedData[:,1]),
            [-8.8, -2.1e-1, -5.0e-2, -5.2e-3, 9.1e-6]
        )
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar, 
            u = fit1[0], 
            a = fit1[1],
            b = fit1[2], 
            c = fit1[3], 
            k = fit1[4]
        ) 
        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-2, 'log'))

    if estar_distribution=="incluster": 
        fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss, 
            np.power(10., ed.inclusterData[:,0]),
            ed.inclusterData[:,1]/max(ed.inclusterData[:,1]),
            [-7.2, -3.9e-1, -1.2e-1, -1.7e-2, 1.6e-4]
            )
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        )
        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-1, 'log'))

    if estar_distribution=="galcenter": 
        fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss, 
            np.power(10., ed.galcenterData[:,0]),
            ed.galcenterData[:,1]/max(ed.galcenterData[:,1]),
            [-4.7, -5.5e-1, -1.8e-1, -1.8e-2, 1.0e-3]
            )
        pdf_unnormed = lambda estar : ed.multigauss(
            y = estar,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        )
        return np.vectorize(ed.generate_invcdf(pdf_unnormed, 1.0e-12, 1.0e-1, 'log'))

fieldData = np.array([[-7.56450, 0.02073], [-7.50032, 0.04945], [-7.45449, 
    0.074], [-7.40865, 0.09855], [-7.36281, 0.12607], [-7.32156, 
    0.15245], [-7.28489, 0.17662], [-7.24363, 0.20363], [-7.20238, 
    0.23064], [-7.16113, 0.25731], [-7.11529, 0.28414], [-7.06945, 
    0.31028], [-7.02361, 0.33597], [-6.97319, 0.36242], [-6.93652, 
    0.38811], [-6.89985, 0.41769], [-6.86318, 0.44728], [-6.82979, 
    0.47512], [-6.76693, 0.50332], [-6.66608, 0.508], [-6.56524, 
    0.5157], [-6.46440, 0.53635], [-6.36356, 0.55731], [-6.26271, 
    0.57744], [-6.15882, 0.58102], [-6.06103, 0.55892], [-5.96935, 
    0.54543], [-5.90518, 0.51636], [-5.85934, 0.48799], [-5.81351, 
    0.45962], [-5.76309, 0.43181], [-5.70808, 0.4068], [-5.65308, 
    0.38199], [-5.59807, 0.35603], [-5.55223, 0.32848], [-5.51556, 
    0.30389], [-5.47889, 0.2793], [-5.44222, 0.25528], [-5.40555, 
    0.23097], [-5.36430, 0.20557], [-5.30929, 0.18175], [-5.24512, 
    0.15767], [-5.17637, 0.13245], [-5.09844, 0.10861], [-5.01135, 
    0.08378], [-4.91509, 0.06145], [-4.81425, 0.04239], [-4.71341, 
    0.0309], [-4.61257, 0.02097], [-4.51172, 0.01984], [-4.33754, 
    0.01724], [-4.24014, 0.01265], [-4.13280, 0.00959], [-4.04051, 
    0.00737], [-3.97542, 0.00759]])

ejectedData = np.array([[-7.47130, 0.00118], [-7.27652, 0.00268], [-7.05391, 
    0.00927], [-6.80348, 0.01868], [-6.56696, 0.04861], [-6.34435, 
    0.05689], [-6.12174, 0.07402], [-5.85739, 0.06913], [-5.63478, 
    0.05671], [-5.39826, 0.04842], [-5.18957, 0.03675], [-4.92522, 
    0.02885], [-4.68870, 0.01736], [-4.46609, 0.01435], [-4.25739, 
    0.01002], [-4.02087, 0.00814], [-3.81217, 0.00344], [-3.57565, 
    0.00231], [-3.32522, 0.00099]])

inclusterData = np.array([[-7.12348, 0.00099], [-6.88696, 0.00249], [-6.59478, 
    0.0008], [-6.27478, 0.0008], [-5.95478, 0.00118], [-5.64870, 
    0.00268], [-5.44000, 0.00551], [-5.20348, 0.01793], [-4.96696, 
    0.03148], [-4.74435, 0.04071], [-4.39652, 0.04071], [-4.29913, 
    0.03694], [-4.03478, 0.02621], [-3.81217, 0.02471], [-3.54783, 
    0.02019], [-3.35304, 0.01096], [-3.13043, 0.01078], [-2.89391, 
    0.00325], [-2.64348, 0.00136]])

galcenterData = np.array([[-4.82673, 0.0186], [-4.52445, 0.16625], [-4.22307, 
    0.21842], [-3.86318, 0.17973], [-3.51015, 0.1252], [-3.17014, 
    0.08535], [-2.84276, 0.05826], [-2.50726, 0.04084], [-2.16264, 
    0.02806], [-1.81864, 0.01836], [-1.50567, 0.01293], [-1.16284, 
    0.00903], [-0.85034, 0.00592], [-0.52809, 0.00435]])

fieldtripleData = np.array([[-7.75, 13], [-7.5, 11], [-7.25, 20], [-7., 
    55], [-6.75, 54], [-6.5, 89], [-6.25, 125], [-6., 
    123], [-5.75, 102], [-5.5, 76], [-5.25, 67], [-5., 
    59], [-4.75, 57], [-4.5, 26], [-4.25, 32], [-4., 14], [-3.75,
     16], [-3.5, 11], [-3.25, 4], [-3., 4], [-2.75, 6], [-2.5, 
    3], [-2.25, 2], [-2., 4], [-1.75, 1], [-1.5, 1], [-1.25, 
    2], [-1., 1], [-0.75, 1], [-0.5, 0], [-0.25, 0], [0., 0]])
