import eccentricdark as ed
import numpy as np
import scipy 

def SnLISA(fp, NModel, AModel): 
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

    return val 

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

def SnLISAdefault(fp): 
    return (ed.SnLISA(fp, "N2", "A5") + ed.SnGal(fp, "N2", "A5"))
