# The "World" class defines an comoving volume filled with binaries
# according to a given mass distribution and number density. 

import os
import numpy as np

class World: 

    
    def __init__(
        self, 
        chi_max, # Units: Mpc
        annual_merger_rate, # Units: binaries * yr^-1 * Mpc^-3  
        buffer_size=10000 # TO DO: Implement this...
    ): 

        self.chi_max = chi_max
        self.annual_merger_rate = annual_merger_rate
    
        self.comoving_volume = (
            (4./3.)
            * np.pi
            * np.power(chi_max, 3.) 
        ) # Units: Mpc^3
    
        print(
            'Initializing world with comoving volume: '
            + f'{self.comoving_volume:.2f}'
            + ' Mpc^3...'
        ) 
            
        self.binary_limit = (
            self.comoving_volume 
            * annual_merger_rate 
        ) # Units: binaries * yr^-1
    
        print(
            'Corresponding to a maximum of '
            + f'{self.binary_limit:.2f}'
            + ' binaries per year.' 
        )
    
        #populate_world(self)
    
    
    def populate_world(
        self,
        mass_distribution,
        mass_args=None
    ):   

        self.mass_distribution = mass_distribution
        self.mass_args = mass_args

        counter = 0
        world_mass = 0. #TO DO: not used... 
        world_full = False
        world_population_m1 = np.zeros(10**9) #TO DO: Improve this 
        world_population_m2 = np.zeros(10**9) #TO DO: Improve this
        world_chi_vals = np.zeros(10**9) #TO DO: Improve this

        while world_full==False: 
            
            if ((counter%int(np.floor(self.binary_limit/100.)))==0): 
                print(
                    'Sample evaluation is '
                    + str(np.floor(100.*counter/self.binary_limit))
                    + '% complete...'
                )
            
            sample = mass_distribution_sampler(
                form=self.mass_distribution, 
                args=mass_args
            )

            sample_m1 = sample[0]
            sample_m2 = sample[1]
            
            M = sample_m1+sample_m2
            
            sample_chi = np.max([10.0**-9, np.random.uniform(0., self.chi_max)])

            if (counter>self.binary_limit):
                world_full=True
            else: 
                world_population_m1[counter] = sample_m1
                world_population_m2[counter] = sample_m2
                world_chi_vals[counter] = sample_chi
                world_mass += M 
                counter += 1

        cut_idx = np.argmin(world_population_m1)

        self.m1_values = world_population_m1[0:cut_idx] #np.array([world_population_m1[0:cut_idx]])
        self.m2_values = world_population_m2[0:cut_idx] #np.array([world_population_m2[0:cut_idx]])
        self.chi_values = world_chi_vals[0:cut_idx] #np.array([world_chi_vals[0:cut_idx]])
        
        max_len=len(self.chi_values)
        print("World contains: ", max_len, " binaries...")
        

    def save(
        self, 
        filepath
    ): 

        if os.path.exists(filepath):
            print("Deleting pre-existing file at path...")
            os.remove(filepath)
        else:
            print("Creating file at path...")
            header_text = (
                'Max comoving distance [Mpc]: ' 
                + f'{self.chi_max:.3f}'
                + '\n'
                + 'Annual Merger Rate [yr^-1 Mpc^-3]: ' 
                + f'{self.annual_merger_rate:.3e}'
                + '\n'
                + 'Mass distribution: ' 
                + self.mass_distribution
                + '\n'
                + 'Mass args: ' 
                + str(self.mass_args)
            )
        
            data = np.stack(
                (self.chi_values, self.m1_values, self.m2_values),
                axis=1
            )

            np.savetxt(filepath, data, header=header_text, comments='#')


def mass_distribution_sampler(
    form, 
    args=None): 

    if form=='flat': 
        return np.random.uniform(
            low=args[0], 
            high=args[1], 
            size=1
        ) #Units M_solar
    elif form=='gaussian': 
        return np.random.normal(
            loc=args[0], 
            scale=args[1], 
            size=1
        ) #Units M_solar      
    elif form=='1602.03842': #Double check this 
        def m1_pdf(m1): #To do: not used... 
            if m1<5: 
                return 0.
            elif m1>100: 
                return 0. 
            elif ((m1>=5) and (m1<=100)): 
                normalization = (
                    1.35 
                    / (np.power(5., -1.35) - np.power(100., -1.35))
                )
                return (normalization * np.power(m1, -2.35))
            else: 
                return False 
            
        def m1_cdf(m1): #To do: not used...
            if m1<5: 
                return 0.
            elif m1>100: 
                return 1. 
            elif ((m1>=5) and (m1<=100)): 
                return (
                    (np.power(5., -1.35) - np.power(m1, -1.35))
                    / (np.power(5., -1.35) - np.power(100., -1.35))
                )
            else: 
                return False 
            
        def m1_invcdf(cdf): 
            return (
                np.power(
                    np.power(5., -1.35) 
                    - cdf * (
                        np.power(5.,-1.35)
                        - np.power(100., -1.35)
                    ), 
                -1./1.35)
            )
        
        random_sample = np.random.random(1)
        m1 = m1_invcdf(random_sample)
        
        m2_max = (100.-m1)
        qmin = (5./m1)
        qmax = 1.

        def q_pdf(q): #To do: not used... 
            if (q<(5./m1)): 
                return 0.
            elif (q>1.):
                return 0
            elif ((q>=(5./m1)) and (q<=1.)): 
                return 1./(qmax-qmin)
            else: 
                return False
        
        def q_cdf(q): #To do: not used...
            if (q<(5./m1)): 
                return 0.
            elif (q>1.):
                return 1
            elif ((q>=(5./m1)) and (q<=1.)): 
                return (1./(qmax-qmin))*(q-qmin)
            else: 
                return False            
            
        def q_invcdf(cdf): 
            return (cdf*(qmax-qmin) + qmin) 
        
        random_sample = np.random.random(1)
        q = q_invcdf(random_sample)
        m2 = m1 * q
        
        return (m1, m2) 
        
    elif form=='1907.02283': 
        def m1_pdf(m1):
            if m1<5: 
                return 0.
            elif m1>50: 
                return 0. 
            elif ((m1>=5) and (m1<=50)): 
                normalization = (
                    1.3 / 
                    (np.power(5., -1.3) - np.power(50., -1.3))
                )
                return (normalization * np.power(m1, -2.3))
            else: 
                return False 
            
        def m1_cdf(m1): 
            if m1<5: 
                return 0.
            elif m1>50: 
                return 1. 
            elif ((m1>=5) and (m1<=50)): 
                return (
                    (np.power(5., -1.3) - np.power(m1, -1.3))
                    / (np.power(5., -1.3) - np.power(50., -1.3))
                )
            else: 
                return False 
            
        def m1_invcdf(cdf): 
            return np.power(
                np.power(5., -1.3) 
                - cdf * (
                    np.power(5.,-1.3)
                    - np.power(50., -1.3)
                ), -1./1.3
            )
        
        random_sample = np.random.random(1)
        m1  = m1_invcdf(random_sample)
        m2 = float(m1)  

        return (m1, m2) 
    
