import os
import scipy
import numpy as np
import eccentricdark as ed

class Cosmology:

    def __init__(
        self,  
        H0=ed.H0, #Units: km/s/Mpcs
        omega_r=ed.omega_r, #Unitless
        omega_m=ed.omega_m, #Unitless
        omega_k=ed.omega_k, #Unitless
        omega_lambda=ed.omega_lambda #Unitless
    ): 

        self.H0 = H0
        self.omega_r = omega_r
        self.omega_m = omega_m
        self.omega_k = omega_k
        self.omega_lambda = omega_lambda
    
    def H(self, z): 
        return ed.H(
            z, 
            self.H0, 
            self.omega_r, 
            self.omega_m, 
            self.omega_k, 
            self.omega_lambda
        )

    def comoving_distance(self, z_i, z_f):
        return ed.comoving_distance(
            z_i, 
            z_f, 
            self.H0, 
            self.omega_r, 
            self.omega_m, 
            self.omega_k, 
            self.omega_lambda
        ) 
    
    def comoving_volume(self, z_i, z_f): 
        return ed.comoving_volume(
            z_i,
            z_f,
            self.H0,
            self.omega_r,
            self.omega_m,
            self.omega_k,
            self.omega_lambda
        )

    def redshift_solver(
        self,
        chi_f, # Units: Mpc
        z_guess # Unitless
    ): 
        return ed.redshift_solver(
            chi_f, 
            z_guess, 
            self.H0,
            self.omega_r,
            self.omega_m,
            self.omega_k,
            self.omega_lambda
        )

    def save(self, filepath): 
        if os.path.exists(filepath):
            print("Deleting pre-existing file at path...")
            os.remove(filepath)
        print("Creating file at path...")
        header_text = (
            'Data columns: H0, omega_r, omega_m, omega_k, omega_lambda'
        )

        data = np.array([
            self.H0,        
            self.omega_r, 
            self.omega_m, 
            self.omega_k, 
            self.omega_lambda
        ])

        np.savetxt(filepath, data.reshape(1, data.shape[0]), header=header_text, comments='#')

def load_cosmology(filepath): 
    if os.path.exists(filepath):
        print("Loading cosmology file at path: "+filepath)
        data = np.loadtxt(filepath)
        H0 = data[0]
        omega_r = data[1]
        omega_m = data[2]
        omega_k = data[3]
        omega_lambda = data[4]
        return ed.Cosmology(H0, omega_r, omega_m, omega_k, omega_lambda)
    else: 
        print("Invalid path, specified file doesn't exist...")

def H(
    z, 
    H0=ed.H0,  #Units: km/s/Mpcs
    omega_r=ed.omega_r, #Unitless
    omega_m=ed.omega_m, #Unitless
    omega_k=ed.omega_k, #Unitless
    omega_lambda=ed.omega_lambda #Unitless
    ): 

    return 100.*np.sqrt(
        omega_r * np.power(1.+z, 4.) 
        + omega_m * np.power(1.+z, 3.) 
        + omega_k * np.power(1.+z, 2.)
        + omega_lambda #Units km/s/Mpc
    )


def comoving_distance(
    z_i, # Unitless
    z_f, # Unitless
    H0=ed.H0,  #Units: km/s/Mpcs
    omega_r=ed.omega_r, #Unitless
    omega_m=ed.omega_m, #Unitless
    omega_k=ed.omega_k, #Unitless
    omega_lambda=ed.omega_lambda #Unitless
    ):

    integrand = lambda z: (
        (ed.c/1000.)
        /H(z, H0, omega_r, omega_m, omega_k, omega_lambda) # Units: Mpc
    )

    return scipy.integrate.quad(integrand, z_i, z_f)[0] # Units: Mpc 


def comoving_volume(
    z_i, # Unitless
    z_f, # Unitless
    H0=ed.H0,  #Units: km/s/Mpcs
    omega_r=ed.omega_r, #Unitless
    omega_m=ed.omega_m, #Unitless
    omega_k=ed.omega_k, #Unitless
    omega_lambda=ed.omega_lambda #Unitless
    ):

    return (
        (4./3.) 
        * np.pi 
        * (
            np.power(
                comoving_distance(
                    0., 
                    z_f, 
                    H0, 
                    omega_r, 
                    omega_m, 
                    omega_k, 
                    omega_lambda), 
                3.
            ) # Units: Mpc^3
            - np.power(
                comoving_distance(
                    0., 
                    z_i, 
                    H0, 
                    omega_r, 
                    omega_m, 
                    omega_k, 
                    omega_lambda), 
                3.
            ) # Units: Mpc^3
        )
    ) #Units: Mpc^3


def redshift_solver(
    chi_f, # Units: Mpc
    z_guess, # Unitless
    H0=ed.H0,  #Units: km/s/Mpcs
    omega_r=ed.omega_r, #Unitless
    omega_m=ed.omega_m, #Unitless
    omega_k=ed.omega_k, #Unitless
    omega_lambda=ed.omega_lambda #Unitless
    ):

    chi_eval = (
        lambda z: (
            chi_f 
            - comoving_distance(
                0, 
                z, 
                H0, 
                omega_r, 
                omega_m, 
                omega_k, 
                omega_lambda)
        )
    )

    return scipy.optimize.fsolve(chi_eval, z_guess) # Units: Mpc



