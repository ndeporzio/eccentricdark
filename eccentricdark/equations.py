import eccentricdark as ed
import numpy as np
import scipy
from scipy import interpolate
import matplotlib.pyplot as plt


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
        * np.power(ed.G_script(e)*(1.-np.power(e, 2.)), -3./2.)
    ) # Unitless


def e_to_fp_interpolator(
    fp_star, 
    e_star,
    e_offset = ed.e_interp_offset_default
    ): 

    bin_count = ed.e_bin_count_default

    e_interp_table = np.geomspace(0.+e_offset, 1.-e_offset, bin_count)
    H_interp_table = H_script(e_interp_table)
    fp_interp_table = (H_script(e_interp_table)/H_script(e_star))*fp_star

    fp_max = fp_interp_table[0]
    fp_min = fp_interp_table[-1]

    fp_of_e_interp = interpolate.interp1d(e_interp_table, fp_interp_table)
    e_of_fp_interp = interpolate.interp1d(fp_interp_table, e_interp_table, bounds_error=False, fill_value=(1., 0.)) 

    return fp_min, fp_max, fp_of_e_interp, e_of_fp_interp 


def e_solver(
    f_p,
    fp_star,
    e_star,
    ): # 1907.02283, eq 4

    #NOTE: Typically better to use e_to_fp_interpolator
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
    if (e==0.): 
        val = (
            np.sqrt(ed.G*m)
            * np.power(a, -3./2.)
            * (1./np.pi) 
        )
    else: 
        val = (
            np.sqrt(ed.G * m)
            * np.power(1.0 + e, ed.gamma)
            * (1.0/np.pi)
            * np.power(a * (1.0 - e**2), -3.0/2.0)
        ) #CAUTION: Blows up for e=1
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
    ) #CAUTION: Blows up for e=1
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
    ) #CAUTION: Blows up for e=1
    return val

# Unit checked 
def tmerge(
    fp, # Units: s^{-1} 
    e, # Units: none
    m1, # Units: kg
    m2  # Units: kg 
):
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
    ) #CAUTION: Blows up for e=1
    
    return val #Units: s

def afe(fp, e, m1, m2): # semi-major axis "a" as a function of "e, fp, m1, m2"
    if e==1.: 
        print("ERROR: afe function blows up for e=1!") 

    m = m1 + m2

    val = (
        (np.power(1.+e, 2.*ed.gamma/3.)*np.power(ed.G * m, 1./3.))
        / ((1.-e**2)*np.power(fp*np.pi, 2./3.)) 
    ) 
    return val 


def BBHevolve(                          
    fp0,                            # Units: s^{-1} (INITIAL frequency of 
                                    # binary at t=t0)
    e0,                             # Units: none (INITIAL eccentricity of 
                                    # binary at t=t0) 
    m1,                             # Units: kg (mass of binary 1) 
    m2,                             # Units: kg (mass of binary 2) 
    t0=(0.*ed.year_in_seconds),     # Units: s (INITIAL time of binary 
                                    # system evolution) 
    tf=(10.*ed.year_in_seconds),    # Units: s (FINAL time of binary system evolution to evolve to) 
    dt=ed.dtevolve,                 # Units: s (time step to take in binary evolution) 
    evolve_factor = 0.001,          # Units: none (upper limit on relative change in 
                                    # a(t) or fp(t) in a time step before time step 
                                    # needs to adaptively be made smaller)  
    verbose=0                       # Units: none (output text verbosity) 
):

    # Calculates black hole binary semimajor axis a(t), eccentricity e(t), and peak emission frequency fp(t)
    # as a function of time given initial state of the system. Time step adaptively becomes smaller if 
    # observable quantities begin to change too quickly. 
    #
    # Output: t, a, e, merger time, merger index, fp 
    # General speed: Fast  

    try: 
        tmerge = ed.tmerge(fp0, e0, m1, m2)
    except:
        print('tmerge blew up')
        tmerge = 1.0*ed.year_in_seconds

    dt = min(dt, tmerge/1.0e3, 1.0*ed.year_in_seconds, tf/1.0e3) 

    m = m1 + m2
    merger = None
    merger_idx = None

    t = [] 
    a = [] 
    e = []  
    fp = [] 
    ap = []
    ep = []

    t.append(t0)
    a.append(ed.afe(fp0, e0, m1, m2))
    e.append(e0)
    fp.append(fp0)
    ap.append(0)
    ep.append(0)

    N = 1
    Rs = (6.*ed.G*m/(ed.c**2)) 
    #Rs = (2.*ed.G*m/(ed.c**2)) 
    while ((a[N-1] > Rs) and (t[-1]< (tf+2.0*dt))):  
        t.append(0)
        a.append(0)
        e.append(0)
        fp.append(0)
        ap.append(0)  
        ep.append(0)  

        good_evolve=False
        while good_evolve==False: 
            t[N] = t[N-1]+dt
            ap[N-1] = ed.dadt(m1, m2, a[N-1], e[N-1])
            ep[N-1] = ed.dedt(m1, m2, a[N-1], e[N-1])             
            a[N] = max(a[N-1] + ap[N-1]*dt, 0.)
            e[N] = max(e[N-1] + ep[N-1]*dt, 0.)
            fp[N] = ed.fpeak(m1, m2, a[N], e[N])

            converge_test = ((a[N-1] - a[N])/a[N-1]) #this ordering to keep quantity positive
            converge_test2 = ((fp[N] - fp[N-1])/fp[N-1]) 

            if (converge_test>evolve_factor): #too large of a change in a[t]
                if (verbose>0): 
                    print("\t Fail: ", converge_test)
                dt = ed.dt_rescale*dt
            elif (converge_test2>evolve_factor): #too large of a change in fp[t]
                if (verbose>0):
                    print("\t Fail: ", converge_test)
                dt = ed.dt_rescale*dt
            else: #good step, continue to next index
                if (verbose>0): 
                    print("a = ", a[N], ", Rs = ", Rs)
                good_evolve=True
                N = len(t)

    if (a[-1]<Rs): 
        merger = t[-1]
        merger_idx = (len(t)-1)
    
    return [
        np.array(t).flatten(), 
        np.array(a).flatten(), 
        np.array(e).flatten(), 
        merger, 
        merger_idx, 
        np.array(fp).flatten()
    ]    

#def fpr(
#    fp0, 
#    e0, 
#    m1, # in kg
#    m2, # in kg
#    ainterp=None, 
#    einterp=None
#): 
#
#
#    if e0==1.:
#        print("ERROR in fpr function. Blows up for e=1!")
# 
#    m = m1 + m2
#
#    if ((ainterp==None) or (einterp==None)): 
#        BBHsol = ed.BBHevolve(fp0, e0, m1, m2)
#        t = BBHsol[0]
#        a = BBHsol[1]
#        e = BBHsol[2]
#        merger = BBHsol[3]
#    
#        ainterp = scipy.interpolate.interp1d(t, a)
#        einterp = scipy.interpolate.interp1d(t, e)
#
#    val = (
#        lambda t : (
#            np.sqrt(ed.G * m)
#            * np.power(1. + einterp(t), ed.gamma)
#            * (1./np.pi)
#            * np.power(ainterp(t) * (1. - np.power(einterp(t), 2.)), -3./2.)  
#        )   
#    ) 
#    
#    return val 

def integrand(
    fp0, # Units: s^{-1} (INITIAL frequency of binary at t=t0) 
    e0, # Units: none (INITIAL eccentricity of binary at t=t0)
    m1, # Units: kg (mass of binary 1) 
    m2, # Units: kg (mass of binary 2) 
    tf, # Units: s (time after entering observation window you wish to define domain of integrand function over) 
    ainterp=None, # Units: m (a(t) interp function, calculated from BBHevolve if not provided)
    einterp=None, # Units: none (e(t) interp function, calculated from BBHevolve if not provided)
    finterp=None, # Units: s^{-1} (fp(t) interp function, calculated from BBHevolve if not provided)
    experiment="LISA"
): 

    # Calculates integrand of SNR integral as a function of time, t, given
    # the initial state of black hole binary system.
    #
    # Output Units: s^{-7/3}   
    # General speed: Fast 

    if ((ainterp==None) or (einterp==None) or (finterp==None)):
        BBHsol = ed.BBHevolve(
            fp0, 
            e0, 
            m1, 
            m2, 
            t0=0.*ed.year_in_seconds,
            tf=tf
        )
        t = BBHsol[0]
        a = BBHsol[1]
        e = BBHsol[2]
        merger = BBHsol[3]
        merger_idx = BBHsol[4]
        f = BBHsol[5]    

        ainterp = scipy.interpolate.interp1d(t, a)
        einterp = scipy.interpolate.interp1d(t, e)
        finterp = scipy.interpolate.interp1d(t, f)
    
    #fpr = ed.fpr(fp0, e0, m1, m2, ainterp, einterp)

    if experiment=="LISA":
        val = (
            lambda t : (
                #np.power(np.pi * fpr(t), 4./3.)
                np.power(np.pi * finterp(t), 4./3.) # Units: s^{-4/3}
                * np.power(1. - einterp(t), 3./2.) # Units: none 
                #*(1./2.32024693e-41) #LISA hn,min 
                #*(1./6.7963438e-49) #DECIGO hn,min
                #* (1./ed.SnLISAdefault(fpr(t)))
                * (1./ed.SnLISAdefault(finterp(t))) # Units: s^{-1} 
            )
        ) # Units: s^{-7/3} 
    elif experiment=="DECIGO":
        try: 
            val = (
                lambda t : (
                    #np.power(np.pi * fpr(t), 4./3.)
                    np.power(np.pi * finterp(t), 4./3.)
                    * np.power(1. - einterp(t), 3./2.)
                    #*(1./2.32024693e-41) #LISA hn,min 
                    #*(1./6.7963438e-49) #DECIGO hn,min
                    #* (1./ed.SnDECIGOdefault(fpr(t)))
                    * (1./ed.SnDECIGOdefault(finterp(t)))
                )
            )
        except:
            pass  
            #print(fpr(t))

    return val # Units: s^{-7/3} (interpolating function for SNR integrand as a function of time) 

def snrintsol(
    fp0, # Units: s^{-1} (INITIAL frequency of binary in detector at t=t0) 
    e0, # Units: none (INITIAL eccentricity of binary in detector at t=t0) 
    m1, # Units: kg (mass of binary 1)
    m2, # Units: kg (mass of binary 2)
    tobserve, # Units: s (max time in seconds that binary is observed after signal first appears in detector)  
    experiment="LISA",
    diagnose=False
): 
     
    BBHsol = ed.BBHevolve(fp0, e0, m1, m2, tf=tobserve)
    t = BBHsol[0]
    a = BBHsol[1]
    e = BBHsol[2]
    tmerger = BBHsol[3]
    tmerger_idx = BBHsol[4]
    f = BBHsol[5]

    if (len(t)<2): 
        return 0. 
    else: 
        ainterp = scipy.interpolate.interp1d(t, a)
        einterp = scipy.interpolate.interp1d(t, e)
        finterp = scipy.interpolate.interp1d(t, f) 
    
        if (tmerger==None): 
            tf = tobserve
        else: 
            tf = min(tobserve, tmerger) 
    
        if (diagnose==True): 
            print("BBHevolve Tmin: ", np.amin(t)/ed.year_in_seconds, " years")
            print("BBHevolve Tmax: ", np.amax(t)/ed.year_in_seconds, " years")
            print("User specified max observe time (ttoday): ", tobserve/ed.year_in_seconds, " years")
            if (tmerger==None): 
                print("No merger happened.")
            else: 
                print("A merger happened, at tmerge = ", tmerger/ed.year_in_seconds, " years")
            print("Integrating BBH evolution with this many time steps: ", len(t))
            print("Between t0 = ", 0., " years and tf = ", tf/ed.year_in_seconds, " years")

        dXdt = []
        X = [0.]
        idx = 0

        integrand = ed.integrand(fp0, e0, m1, m2, tf, ainterp, einterp, finterp, experiment)
        while t[idx]<tf:
            if (diagnose==True): 
                print("Integrating time ", t[idx]/ed.year_in_seconds, " of ", tf/ed.year_in_seconds, " years. X = ", X[-1]) 
            dXdt.append(integrand(t[idx]))
            X.append(X[idx] + dXdt[idx]*(t[idx+1]-t[idx])) # Units: s^{-4/3} 
            idx+=1
     
        return X[-1] # Units: s^{-4/3} 


def SNR_LISA(
    r, # Units: m (distance of binary system from the detector) 
    fp0, # Units: s^{-1} (INITIAL frequency of binary in detector at t=t0) 
    e0, # Units: none (INITIAL eccentricity of binary in detector at t=t0)
    m1, # Units: kg (mass of binary 1)
    m2, # Units: kg (mass of binary 2)
    tobserve, # Units: s (max time in seconds that binary is observed after signal first appears in detector)   
    experiment="LISA"
): 
    
    m = m1 + m2 # Units: kg
    mu = ed.m_reduced(m1, m2) # Units: kg  

    val = np.sqrt(
        (2.*64./5.) 
        * np.power(ed.c, -8.) 
        * np.power(r, -2.) 
        * np.power(ed.G * np.power(mu, 3./5.) * np.power(m, 2./5.), 10./3.)
        * ed.snrintsol(fp0, e0, m1, m2, tobserve, experiment)
    )

    return val # Units: none 


def roffmSNR8(
    r, # Units: Mpc (distance of binary system from the detector)  
    fp0, # Units: s^{-1} (INITIAL frequency of binary in detector at t=t0)
    e0,  # Units: none (INITIAL eccentricity of binary in detector at t=t0)
    m1,  # Units: kg (mass of binary 1)
    m2, # Units: kg (mass of binary 2)
    tobserve, # Units: years (max time in seconds that binary is observed after signal first appears in detector)   
    experiment="LISA"
    ): 

    # Return the comoving distance at which a binary system signal
    # with specified masses, signal eccentricity, signal peak 
    # frequency and observation time can be measured with SNR>8 

    #print('SNR', ed.SNR_LISA(r*(10.**6)*ed.pc_in_meters, fp, e, m1, m2, t*ed.year_in_seconds, experiment))
    #print('Params')
    #print(r*(10.**6)*ed.pc_in_meters, fp, e, m1, m2, t*ed.year_in_seconds)

    SNRVal = 8.

    return (
        # This only works because SNR scales linearly with distance 
        # i.e. deccreasing distance by factor of 2 increases SNR by a factor of 2 
        (r * ((10.**6) * ed.pc_in_meters)) # Units: m (distance of source in meters)       
        * (1./SNRVal) # Units: none 
        * ed.SNR_LISA(
            r*(10.**6)*ed.pc_in_meters, 
            fp0, 
            e0, 
            m1,
            m2, 
            tobserve*ed.year_in_seconds, 
            experiment
        ) # Units: none 
    ) # Units: m (distance at which observation of system can happen with SNR=8) 

