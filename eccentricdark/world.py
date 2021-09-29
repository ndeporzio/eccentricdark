# The "World" class defines an comoving volume filled with binaries
# according to a given mass distribution and number density. 

import os
import numpy as np
import eccentricdark as ed

class World: 

    
    def __init__(
        self, 
        chi_max, # Units: Mpc
        annual_merger_rate, # Units: binaries * yr^-1 * Mpc^-3  
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
    
    
    def populate_world(
        self,
        mass_distribution,
        mass_args=None
    ):   

        self.mass_distribution = mass_distribution
        self.mass_args = mass_args

        min_chi = 1.0e-9

        counter = 0
        world_full = False
        world_mass = 0.

        world_population_m1 = np.zeros(np.int(np.floor(self.binary_limit)))  
        world_population_m2 = np.zeros(np.int(np.floor(self.binary_limit))) 
        world_chi_vals = np.zeros(np.int(np.floor(self.binary_limit))) 
        world_mc_vals = np.zeros(np.int(np.floor(self.binary_limit)))

        while world_full==False: 
            
            if ((counter%int(np.floor(self.binary_limit/100.)))==0): 
                print(
                    'Sample evaluation is '
                    + str(np.floor(100.*counter/self.binary_limit))
                    + '% complete...'
                )
            
            sample = ed.mass_distribution_sampler(
                form=self.mass_distribution, 
                args=mass_args
            )

            sample_m1 = sample[0]
            sample_m2 = sample[1]
            
            M = sample_m1+sample_m2
            mc = ed.m_chirp(sample_m1, sample_m2) 
            
            sample_chi = np.random.uniform(min_chi, self.chi_max)

            if (counter>=np.int(np.floor(self.binary_limit))):
                world_full=True
            else: 
                world_population_m1[counter] = sample_m1
                world_population_m2[counter] = sample_m2
                world_chi_vals[counter] = sample_chi
                world_mc_vals[counter] = mc
                world_mass += M 
                counter += 1

        cut_idx = np.argmin(world_population_m1)

        self.world_mass = world_mass
        self.m1_values = world_population_m1
        self.m2_values = world_population_m2
        self.chi_values = world_chi_vals
        self.mc_values = world_mc_vals
        
        max_len=len(self.chi_values)
        print("World contains: ", max_len, " binaries...")
        

    def save(
        self, 
        filepath
    ): 

        if os.path.exists(filepath):
            print("Deleting pre-existing file at path...")
            os.remove(filepath)

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
            + '\n'
            + 'Data columns: chi, m1, m2, mc'
        )
        
        data = np.stack(
            (self.chi_values, self.m1_values, self.m2_values, self.mc_values),
            axis=1
        )

        np.savetxt(filepath, data, header=header_text, comments='#')

def load_world(filepath): 
    header_text = []

    if os.path.exists(filepath):
        print("Loading world file at path: "+filepath)     
        with open(filepath, 'r') as f: 
            lines = f.readlines()
        for line in lines: 
            if (line[0]=='#'): 
                header_text.append(line)

        chi_max = float(header_text[0].split(
            'Max comoving distance [Mpc]: ')[1])
        annual_merger_rate = float(header_text[1].split(
            'Annual Merger Rate [yr^-1 Mpc^-3]: ')[1])
        mass_distribution = header_text[2].split(
            'Mass distribution: ')[1]
        mass_args = header_text[3].split(
            'Mass args: ')[1]
        if (mass_args[0:4]=='None'): 
            mass_args=None
        else: 
            mass_args = np.array(
                mass_args[1:-1].split(', '), 
                dtype='float'
            )

        data = np.loadtxt(filepath)

        w = ed.World(chi_max, annual_merger_rate)
        w.mass_distribution = mass_distribution
        w.mass_args = mass_args
        w.chi_values = data[:, 0]
        w.m1_values = data[:, 1]
        w.m2_values = data[:, 2]
        w.mc_values = data[:, 3]
        w.world_mass = np.sum(w.m1_values + w.m2_values)

        print("World loaded: ")
        print("\t Max comoving distance [Mpc]: ", chi_max)
        print("\t Annual Merger Rate [yr^-1 Mpc^-3]: ", annual_merger_rate)
        print("\t Mass distribution: ", mass_distribution)
        print("\t Mass args: ", mass_args)

        return w
      
    else: 
        print("Invalid filepath...")
        

