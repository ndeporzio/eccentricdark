import eccentricdark as ed
import numpy as np
import scipy
from scipy import interpolate

def m_chirp(m1, m2): 
    if ((m1==0) or (m2==0)): 
        return 0.
    else: 
        return np.power(m1*m2, 3./5.)/np.power(m1+m2, 1./5.) # Units: mass^1

def m_reduced(m1, m2): 
    if ((m1==0.) or (m2==0.)): 
        return 0.
    else: 
        return ((m1*m2)/(m1+m2))

def G_script( # 1907.02283, eq 2
    e #Unitless
    ): 
    return (
        (np.power(e, 12./19.)/(1-np.power(e, 2.))) 
        * np.power(1. + (121./304.)*np.power(e, 2.), 870./2299.) 
    ) #Unitless

def H_script( # 1907.02283, eq 4
    e #Unitless
    ):
    return (
        np.power(1.+e, ed.gamma) 
        / np.power((1.-np.power(e, 2.))*G_script(e), 3./2.)
    ) # Unitless


def e_to_fp_interpolator(
    fp_star, 
    e_star,
    e_offset = ed.e_interp_offset_default
    ): 

    bin_count = ed.e_bin_count_default

    e_interp_table = np.linspace(0.+e_offset, 1.-e_offset, bin_count)
    H_interp_table = H_script(e_interp_table)
    fp_interp_table = (H_script(e_interp_table)/H_script(e_star))*fp_star

    fp_max = fp_interp_table[0]
    fp_min = fp_interp_table[-1]

    fp_of_e_interp = interpolate.interp1d(e_interp_table, fp_interp_table)
    e_of_fp_interp = interpolate.interp1d(fp_interp_table, e_interp_table) 

    return fp_min, fp_max, fp_of_e_interp, e_of_fp_interp 

def e_solver(
    f_p,
    fp_star,
    e_star,
    ): # 1907.02283, eq 4

    e_offset = ed.e_interp_offset_default
    e_bin_count = ed.e_bin_count_default

    e_interp_table = np.linspace(0.+e_offset, 1.-e_offset, e_bin_count)
    H_interp_table = H_script(e_interp_table)
    H_interp_func = interpolate.interp1d(e_interp_table, H_interp_table)
    H_inv_interp_table = H_interp_table - ((f_p/fp_star) * H_script(e_star)) 
    H_inv_interp_func = interpolate.interp1d(e_interp_table, H_inv_interp_table)

    if H_inv_interp_table[-1]>0:  
        print('This peak frequency never occurs - eccentricity saturated at e=1!')
        return 1
    elif H_inv_interp_table[0]<0: 
        print('This peak frequency never occurs - eccentricity saturated at e=0!')
        return 0
    else: 
        solverfunc = lambda x: H_inv_interp_func(x)
        e_val = scipy.optimize.bisect(solverfunc, 0.+e_offset, 1.-e_offset)
        return e_val # Unitless

def e_solver_2(
    fp,
    fp_max, 
    fp_min, 
    e_of_fp):

    if (fp < fp_min): 
        print('This peak frequency never occurs - eccentricity saturated at e=1!')
        return 1
    elif (fp > fp_max): 
        print('This peak frequency never occurs - eccentricity saturated at e=0!')
        return 0
    else: 
        return e_of_fp(fp) 


def F_script(e): # 1907.02283, eq 6 
    return (
        (
            np.power(1.+e, (8.*ed.gamma/3.)-(1./2.))
            / np.power(1.-e, 3./2.)
        )
        * np.power(
            ((1.+e)*(1.+(7./8.)*(e**2.)))
            - ((ed.gamma/288.)*e*(304.+121.*(e**2.)))
        , -1.)
    )

def dtdfp(
    mc, # in Msun  
    fp, 
    e): # 1907.02283, eq 5
    return (
        (
            (5.*np.power(ed.c, 5.)) 
            / (96*np.power(np.pi, 8./3.))
        )
        * np.power(ed.G * mc * ed.msun_in_kg, -5./3.)
        * np.power(fp, -11./3.)
        * ed.F_script(e) 
    )

def dNdfp_estarfixed_mcfixed(R, chimax, mc, fp, fpstar, estar): 
    e = ed.e_solver(fp, fpstar, estar)

    return (
        R
        * ((4./3.)*np.pi*np.power(chimax, 3.))
        * ed.dtdfp(mc, fp, e)
    )

def N_estarfixed(R, f, chi_max, mass_distribution, estar, fpstar): 
    integrand  = lambda mc, fp, r: (
        R
        * ed.mc_pdf(mc, mass_distribution)
        * (4.0 * np.pi * np.power(r, 2.))
        * ed.dtdfp(mc, fp, e_solver(fp, fpstar, estar))
    ) 

    return scipy.integrate.tplquad(
        integrand, 
        np.power(2.0, -0.2)*5.0*ed.msun_in_kg, 
        np.power(2.0, -0.2)*50.0*ed.msun_in_kg,
        lambda x: f, 
        lambda x: f*np.power(10.0, 0.1), 
        lambda x, y: 0., 
        lambda x, y: chi_max*np.power(10., 6.)*ed.pc_in_meters
    )


def fpeak(m1, m2, a, e): # equation 3
    m = m1 + m2
    val = (
        np.sqrt(ed.G * m)
        * np.power(1.0 + e, ed.gamma)
        * (1.0/np.pi)
        * np.power(a * (1.0 - e**2), -3.0/2.0)
    )
    return val 

def dadt(m1, m2, a, e): # equation 1a 
    m = m1 + m2
    mu = ed.m_reduced(m1, m2) 
    
    val = (
        -1.0
        * (64./5.)
        * np.power(ed.G, 3.)
        * mu
        * np.power(m, 2.)
        * np.power(ed.c, -5.)
        * np.power(a, -3.)
        * (1.0 + (73./24.)*(e**2) + (37./96.)*(e**4))
        * np.power(1.0 - e**2, -7./2.)
    )
    return val 

def dedt(m1, m2, a, e): # equation 1b
    m = m1 + m2
    mu = ed.m_reduced(m1, m2)
    
    val = (
        -1.0
        * (304./15.)
        * np.power(ed.G, 3.)
        * mu
        * np.power(m, 2.)
        * np.power(ed.c, -5.)
        * np.power(a, -4.)
        * e
        * (1.0 + (121./304.)*(e**2))
        * np.power(1.0 - e**2, -5./2.)
    )
    return val

def tmerge(fp, e, m1, m2):
    m = m1 + m2
    mu = ed.m_reduced(m1, m2) 

    val = (
        (5./256.)
        * np.power(ed.c, 5.)
        * np.power(ed.G, -3.)
        * (1./mu)
        * np.power(m, -2.) 
        * np.power(np.power(1.+e, 2.*ed.gamma/3.)*np.power(ed.G * m, 1./3.) , 4.)
        * np.power((1.-e**2)*np.power(fp*np.pi, 2./3.), -4.) 
        * np.power(1.0 - e**2 , 7./2.)
    )
    
    return val 

def afe(fp, e, m1, m2): # semi-major axis "a" as a function of "e, fp, m1, m2"
    m = m1 + m2

    val = (
        (np.power(1.+e, 2.*ed.gamma/3.)*np.power(ed.G * m, 1./3.))
        / ((1.-e**2)*np.power(fp*np.pi, 2./3.)) 
    ) 
    return val 


def BBHevolve(                          #Fast 
    fp0, e0, m1, m2, 
    t0=(0.*ed.year_in_seconds), 
    tf=(10.*ed.year_in_seconds), 
    dt=ed.dtevolve
):
    m = m1 + m2
    merger = False
    size = int((tf-t0)/dt)+1

    t = size * [0]
    a = size * [0]
    e = size * [0]
    ap = size * [0]
    ep = size * [0]

    t[0] = t0
    a[0] = ed.afe(fp0, e0, m1, m2)
    e[0] = e0

    t[-1] = tf

    for idx in range(len(t)):
        if (idx!=0) and (idx!=(len(t)-1)):
            t[idx]=t[idx-1]+dt

    for idx in range(len(t)):
        ap[idx] = ed.dadt(m1, m2, a[idx], e[idx])
        ep[idx] = ed.dedt(m1, m2, a[idx], e[idx])

        if (idx!=(len(t)-1)): #do for all but last entry 
            a[idx+1] = a[idx] + ap[idx]*dt
            e[idx+1] = e[idx] + ep[idx]*dt

            if (a[idx+1] < (6.*ed.G*m/(ed.c**2))): 
                merger = t[idx]
                break 
                

    return [np.array(t).flatten(), np.array(a).flatten(), np.array(e).flatten(), merger]    

def fpr(fp0, e0, m1, m2, ainterp=None, einterp=None): 

    m = m1 + m2

    if ((ainterp==None) or (einterp==None)): 
        BBHsol = ed.BBHevolve(fp0, e0, m1, m2)
        t = BBHsol[0]
        a = BBHsol[1]
        e = BBHsol[2]
        merger = BBHsol[3]
    
        ainterp = scipy.interpolate.interp1d(t, a)
        einterp = scipy.interpolate.interp1d(t, e)

    val = (
        lambda t : (
            np.sqrt(ed.G * m)
            * np.power(1. + einterp(t), ed.gamma)
            * (1./np.pi)
            * np.power(ainterp(t) * (1. - np.power(einterp(t), 2.)), -3./2.)  
        )   
    ) 
    
    return val 

def integrand(fp0, e0, m1, m2, tf=(10.*ed.year_in_seconds), ainterp=None, einterp=None): #Fast 

    if ((ainterp==None) or (einterp==None)):
        BBHsol = ed.BBHevolve(fp0, e0, m1, m2, tf=tf)
        t = BBHsol[0]
        a = BBHsol[1]
        e = BBHsol[2]
        merger = BBHsol[3]
    
        ainterp = scipy.interpolate.interp1d(t, a)
        einterp = scipy.interpolate.interp1d(t, e)

    fpr = ed.fpr(fp0, e0, m1, m2, ainterp, einterp)

    val = (
        lambda t : (
            np.power(np.pi * fpr(t), 4./3.)
            * np.power(1. - einterp(t), 3./2.)
            * (1./ed.SnLISAdefault(fpr(t)))
        )
    )

    return val 

def snrintsol(
    fp0, # signal frequency in detector
    e0, # signal eccentricity in detector
    m1, # binary mass 1 in Msun
    m2, # binary mass 2 in Msun 
    ttoday # observation time in seconds 
): #Slow
     
    BBHsol = ed.BBHevolve(fp0, e0, m1, m2, tf=ttoday)
    t = BBHsol[0]
    a = BBHsol[1]
    e = BBHsol[2]
    merger = BBHsol[3]

    ainterp = scipy.interpolate.interp1d(t, a)
    einterp = scipy.interpolate.interp1d(t, e)

    if (merger==False): 
        tf = ttoday 
    else: 
        tf = merger 

    dt = t[1]-t[0]
    t0 = (0. * ed.year_in_seconds)

    size = int((tf-t0)/dt)+1

    t = size * [0]
    sol = size * [0]
    solp = size * [0]

    t[0] = t0
    sol[0] = 0.

    t[-1] = tf

    for idx in range(len(t)):
        if (idx!=0) and (idx!=(len(t)-1)):
            t[idx]=t[idx-1]+dt

        solp[idx] = ed.integrand(fp0, e0, m1, m2, tf, ainterp, einterp)(t[idx])

        if (idx!=(len(t)-1)): #do for all but last entry 
            sol[idx+1] = sol[idx] + solp[idx]*dt

    return sol[-1] 
    

def SNR_LISA(
    r, # distance in pc 
    fp0, # frequency in detector 
    e0, # eccentricity in detector
    m1, # binary mass 1 in Msun
    m2, # binary mass 2 in Msun
    ttoday # observation time in seconds
): 
    
    m = m1 + m2
    mu = ed.m_reduced(m1, m2) 

    val = np.sqrt(
        (2.*64./5.) 
        * np.power(ed.c, -8.) 
        * np.power(r, -2.) 
        * np.power(ed.G * np.power(mu, 3./5.) * np.power(m, 2./5.), 10./3.)
        * ed.snrintsol(fp0, e0, m1, m2, ttoday)
    )

    return val 


def roffmSNR8(
    r, # distance in Mpc 
    fp, # peak frequency of signal 
    e,  # Eccentricity of signal 
    m1,  # Binary 1 mass in Msun
    m2, # Binary 2 mass in Msun
    t #years
    ): 

    # Return the comoving distance at which a binary system signal
    # with specified masses, signal eccentricity, signal peak 
    # frequency and observation time can be measured with SNR>8 

    return (
        (r/8.)
        * ((10.**6) * ed.pc_in_meters)
        * ed.SNR_LISA(r*(10.**6)*ed.pc_in_meters, fp, e, m1, m2, t*ed.year_in_seconds)
    )

