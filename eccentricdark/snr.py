import eccentricdark as ed
import numpy as np
import scipy 

def SnAcc(fp, NModel): 
    if (NModel=="N1"):
        return (
            9.0e-28
            * np.power(2.0 * np.pi * fp, -4) 
            * (1.0 + 1.0e-4 * np.power(fp, -1))
        )
    elif (NModel=="N2"): 
        return (
            9.0e-30
            * np.power(2.0 * np.pi * fp, -4)
            * (1.0 + 1.0e-4 * np.power(fp, -1))
        )

def SnSn(fp, AModel): 
    if (AModel=="A1"): 
        return 1.98e-23
    if (AModel=="A2"):
        return 2.22e-23
    if (AModel=="A5"):
        return 2.96e-23

def SnOmn(fp): 
    return 2.65e-23 

def SnLISA(
    fp, # Units: s^{-1}
    NModel, # String - specify LISA model
    AModel # String - specify LISA model
): 
    if (AModel=="A1"): 
        L = 1.0e9
    elif (AModel=="A2"): 
        L = 2.0e9
    elif (AModel=="A5"): 
        L = 5.0e9

    val = (
        (20./3.)
        * (
            (4.0 * ed.SnAcc(fp, NModel))
            + ed.SnSn(fp, AModel)
            + ed.SnOmn(fp) 
        )
        * np.power(L, -2.)
        * (1. + np.power(fp/(0.41*ed.c/(2.0*L)), 2.))  
    ) 

    return val # Units: s (should be)

def SnGal(fp, NModel, AModel): 
    if (AModel=="A1"): 
        if ((fp >= 1.0e-5) and (fp < 5.3e-4)):
            return (1.55206e-43 * np.power(fp, -2.1))
        elif ((fp >= 5.3e-4) and (fp < 2.2e-3)): 
            return (2.9714e-47 * np.power(fp, -3.235))
        elif ((fp >= 2.2e-3) and (fp < 4.0e-3)):        
            return (1.517e-51 * np.power(fp, -4.85))
        elif ((fp >= 4.0e-3) and (fp < 5.3e-3)):        
            return (6.706e-58 * np.power(fp, -7.5))
        elif ((fp >= 5.3e-3) and (fp < 1.0e-2)):        
            return (2.39835e-86 * np.power(fp, -20.0))
        else:
            return 0.
    elif (AModel=="A5"):
        if ((fp >= 10.**-5) and (fp < 10.**-3)):        
            return (10**-44.62 * np.power(fp, -2.3))
        elif ((fp >= 10.**-3) and (fp < 10.**-2.7)):        
            return (10.**-50.92 * np.power(fp, -4.4))
        elif ((fp >= 10.**-2.7) and (fp < 10.**-2.4)):        
            return (10.**-62.8 * np.power(fp, -8.8))
        elif ((fp >= 10.**-2.4) and (fp < 10.**-2)):        
            return (10.**-89.68 * np.power(fp, -20.))
        else:
            return 0. 

    # output # Units: s (should be) 

def SnLISAdefault( # Output LISA power spectral density 
    fp # Units: s^{-1}
): 
    return ed.SnLISA(fp, "N2", "A5") + ed.SnGal(fp, "N2", "A5") # Units: s  

def SnDECIGOdefault(fp):
    fp_vals = np.array([2.14350692e-03, 2.43759044e-03, 2.77202145e-03, 3.08550647e-03,
       3.43444317e-03, 3.82284077e-03, 4.25516184e-03, 4.83895888e-03,
       5.38619169e-03, 5.99531049e-03, 6.67331389e-03, 7.42799200e-03,
       8.44709300e-03, 9.40236594e-03, 1.04656697e-02, 1.16492214e-02,
       1.29666196e-02, 1.47456058e-02, 1.64131710e-02, 1.82693194e-02,
       2.03353776e-02, 2.26350843e-02, 2.57405584e-02, 2.86515313e-02,
       3.18917031e-02, 3.54983025e-02, 3.95127685e-02, 4.49338167e-02,
       5.00153351e-02, 5.56715172e-02, 6.19673509e-02, 6.89751740e-02,
       7.84383869e-02, 8.73089022e-02, 9.71825748e-02, 1.10515769e-01,
       1.23013882e-01, 1.36925394e-01, 1.55711198e-01, 1.77074366e-01,
       2.01368504e-01, 2.28995734e-01, 2.60413348e-01, 2.96141376e-01,
       3.36771196e-01, 3.82975320e-01, 4.35518529e-01, 4.95270528e-01,
       5.63220344e-01, 6.40492697e-01, 7.28366614e-01, 8.28296602e-01,
       9.41936722e-01, 1.07116797e+00, 1.21812940e+00, 1.38525356e+00,
       1.57530671e+00, 1.79143466e+00, 2.03721479e+00, 2.31671531e+00,
       2.63456257e+00, 2.99601765e+00, 3.40706340e+00, 3.87450356e+00,
       4.40607528e+00, 5.01057724e+00, 5.69801529e+00, 6.47976803e+00,
       7.05965305e+00])
    root_hn_over_fp = np.array([7.10591460e-22, 5.88461425e-22, 4.92183041e-22, 4.17831557e-22,
       3.54711958e-22, 3.01127502e-22, 2.55637766e-22, 2.13812780e-22,
       1.81513216e-22, 1.54092977e-22, 1.30814967e-22, 1.11053443e-22,
       9.28839498e-23, 7.88524637e-23, 6.69406399e-23, 5.68282722e-23,
       4.82435262e-23, 4.03503857e-23, 3.42548667e-23, 2.90801655e-23,
       2.46871790e-23, 2.09578177e-23, 1.75289017e-23, 1.48809034e-23,
       1.26329242e-23, 1.07245353e-23, 9.10443656e-24, 7.61485644e-24,
       6.46452043e-24, 5.48795959e-24, 4.65892261e-24, 3.95512386e-24,
       3.30802463e-24, 2.85042274e-24, 2.38406405e-24, 1.99400648e-24,
       1.69278249e-24, 1.43706281e-24, 1.20194445e-24, 1.02544956e-24,
       8.92412147e-25, 7.68964005e-25, 6.89428246e-25, 6.36801425e-25,
       5.99984791e-25, 5.88191820e-25, 5.94059043e-25, 6.18119058e-25,
       6.43153529e-25, 6.69201923e-25, 7.10265901e-25, 7.53849671e-25,
       8.00107855e-25, 8.66230718e-25, 9.19384897e-25, 9.75800755e-25,
       1.05644331e-24, 1.12126943e-24, 1.21393386e-24, 1.31425632e-24,
       1.40881669e-24, 1.52524474e-24, 1.66776637e-24, 1.82360548e-24,
       1.99400648e-24, 2.22404466e-24, 2.45612130e-24, 2.76679754e-24,
       2.85042274e-24])
    hn = np.power(root_hn_over_fp, 2.)/fp_vals

    if (fp>fp_vals[-1]): 
        return np.power(fp/fp_vals[-1], 2.)*hn[-1]+ed.SnGal(fp, "N2", "A5")
    elif (fp<fp_vals[0]): 
        return np.power(fp/fp_vals[0], -5.)*hn[0]+ed.SnGal(fp, "N2", "A5")
    else: 
        return (
            np.power(10., scipy.interpolate.interp1d(np.log10(fp_vals), np.log10(hn))(np.log10(fp)))
            + ed.SnGal(fp, "N2", "A5")
        )
   
