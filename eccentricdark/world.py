# The "World" class defines an comoving volume filled with binaries
# according to a given mass distribution and number density. 

import os
import numpy as np
import scipy
import eccentricdark as ed

class World: 

    
    def __init__(
        self, 
        chi_max, # Units: Mpc
        annual_merger_rate, # Units: binaries * yr^-1 * Mpc^-3  
    ): 

        self.world_mass=None
        self.is_evolved=False
        self.use_cosmology=False

        self.estar_distribution=None
        self.estar_args=None
        self.estar=None
        self.estar_invcdf=None

        self.e_of_fp_interp=None

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
        mass_args=None, 
        chi_distribution=None, 
        chi_args=None,
        cosmo=False
    ):   

        self.mass_distribution = mass_distribution
        self.mass_args = mass_args
        self.chi_distribution = chi_distribution
        self.chi_args = chi_args

        if cosmo==True: 
            self.use_cosmology=True

        counter = 0
        world_full = False
        world_mass = 0.

        world_population_m1 = np.zeros(np.int(np.floor(self.binary_limit)))  
        world_population_m2 = np.zeros(np.int(np.floor(self.binary_limit))) 
        world_chi_vals = np.zeros(np.int(np.floor(self.binary_limit))) 
        world_mc_vals = np.zeros(np.int(np.floor(self.binary_limit)))

        while world_full==False: 
            
            #if ((counter%int(np.floor(self.binary_limit/100.)))==0): 
            #    print(
            #        'Sample evaluation is '
            #        + str(np.floor(100.*counter/self.binary_limit))
            #        + '% complete...', end=""
            #    )
            #    print("\r", end="")
            
            sample = ed.mass_distribution_sampler(
                form=self.mass_distribution, 
                args=mass_args
            )

            sample_m1 = sample[0]
            sample_m2 = sample[1]
            
            M = sample_m1+sample_m2
            mc = ed.m_chirp(sample_m1, sample_m2)

            sample_chi = ed.chi_distribution_sampler(
                chi_distribution, 
                chi_args
            )

            if (counter>=np.int(np.floor(self.binary_limit))):
                world_full=True
                print("Sample evaluation is 100% complete.")
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

        if self.use_cosmology==True: 
            self.redshifts = np.zeros(len(self.chi_values))
            cosmo = ed.Cosmology()
            for idx, val in enumerate(self.m1_values):
                self.redshifts[idx] = cosmo.redshift_solver(self.chi_values[idx], 0.1)
                self.mc_values[idx] = (1.+self.redshifts[idx])*self.mc_values[idx]
        
    def initialize_eccentricities(
        self,
        estar_distribution,
        args=None
    ):

        if (self.world_mass==None): 
            print("Need to 'populate_world' masses first...")
            return

        self.estar_distribution = estar_distribution
        self.estar_args = args

        print("Building eccentricity interpolator...")
        self.estar_invcdf = ed.estar_sampler(estar_distribution, args) 
        self.estar_sampler = (lambda : self.estar_invcdf(np.random.random(1))) 

        print("Generating eccentricity samples...") 
        #self.estar = self.estar_invcdf(np.random.random(len(self.mc_values)))
        self.estar = [0]*(len(self.mc_values))

        for mc_idx, mc_val in enumerate(self.mc_values):
            self.estar[mc_idx] = self.estar_sampler()[0]


    def solve_evolution(
        self,
        mode, 
        fpstar=10.
    ): 
        if (self.world_mass==None):
            print("Need to 'populate_world' masses first...")
            return
    
        self.fp_min = np.zeros(len(self.mc_values))
        self.fp_max = np.zeros(len(self.mc_values))
        self.fp_of_e_interp = [0]*len(self.mc_values)
        self.e_of_fp_interp = [0]*len(self.mc_values)

        if (mode=="fixed"): 
            if np.any(self.estar != self.estar[0]): 
                print("Cannot evaluate e[fp] in 'fixed' mode if binaries "
                    + "have different e* values...")
                return
            
            evolution = ed.e_to_fp_interpolator(fp_star=fpstar, e_star=self.estar[0])
            for mc_idx, mc_val in enumerate(self.mc_values):
                self.fp_min[mc_idx] = evolution[0]
                self.fp_max[mc_idx] = evolution[1]
                self.fp_of_e_interp[mc_idx] = evolution[2]
                self.e_of_fp_interp[mc_idx] = evolution[3]

        elif (mode=="individual"): # Warning: high memory use
            if (self.estar_distribution=='fixed'):
                print("The e* distribution is constant. Re-run using 'fixed' mode...")
                return 
                
            print("Evaluating e[fp] for given e*, fp*...")
            for mc_idx, mc_val in enumerate(self.mc_values):
                if (mc_idx%10==0): 
                    print("Solving ", mc_idx+1, "/", len(self.mc_values), " ...", end="")
                    print("\r", end="") 
                evolution = ed.e_to_fp_interpolator(fp_star=fpstar, e_star=self.estar[mc_idx])
                self.fp_min[mc_idx] = evolution[0]
                self.fp_max[mc_idx] = evolution[1]
                self.fp_of_e_interp[mc_idx] = evolution[2]
                self.e_of_fp_interp[mc_idx] = evolution[3]
        elif (mode=="interpolate"): 
            if (self.estar_distribution=='fixed'):
                print("The e* distribution is constant. Re-run using 'fixed' mode...")
                return
            if (self.is_evolved==False): 
                estar_max = np.max(self.estar)
                estar_min = np.min(self.estar)
                de = 0.01
    
                estar_table = np.geomspace((1.-de)*estar_min, (1.+de)*estar_max, 100)
    
                fp_min_table = np.zeros(len(estar_table))
                fp_max_table = np.zeros(len(estar_table))
                fp_of_e_interp_table = [0]*len(estar_table)
                e_of_fp_interp_table = [0]*len(estar_table)
            
                def e_evaluator(e_idx, e_of_fp, fp):
                    if (fp < fp_min_table[e_idx]): 
                        return 1.
                    elif (fp > fp_max_table[e_idx]): 
                        return 0.
                    else: 
                         return e_of_fp(fp) 
    
                for e_idx, e_val in enumerate(estar_table):
                    if (e_idx%10==0): 
                        print("Building interpolation ", e_idx, "/", len(estar_table), "...", end="")
                        print("\r", end="")
                    evolution = ed.e_to_fp_interpolator(fp_star=10., e_star=e_val)
                    fp_min_table[e_idx] = evolution[0]
                    fp_max_table[e_idx] = evolution[1]
                    fp_of_e_interp_table[e_idx] = evolution[2](np.linspace(1.e-9, 1.-1e-9, 100))
                    e_of_fp_interp_table[e_idx] = np.array([e_evaluator(
                        e_idx,
                        evolution[3],
                        fp
                    ) for fp in np.logspace(-6, 2, 801)])
            
                self.fp_min_of_estar_interp = scipy.interpolate.interp1d(estar_table, fp_min_table)
                self.fp_max_of_estar_interp = scipy.interpolate.interp1d(estar_table, fp_max_table)
                self.fp_of_e_for_estar_interp = scipy.interpolate.interp2d(
                    np.linspace(1.e-9, 1.-1e-9, 100),
                    estar_table,
                    fp_of_e_interp_table
                )
                self.e_of_fp_for_estar_interp = scipy.interpolate.interp2d(
                    np.logspace(-6, 2, 801),
                    estar_table, 
                    e_of_fp_interp_table
                ) 

                self.is_evolved=True

            print("Evaluating e[fp] for given e*, fp*...")
            for mc_idx, mc_val in enumerate(self.mc_values):
                if (mc_idx%10==0):
                    print("Solving ", mc_idx, "/", len(self.mc_values)-1, " ...", end="")
                    print("\r", end="")
                    
                self.fp_min[mc_idx] = self.fp_min_of_estar_interp(self.estar[mc_idx])
                self.fp_max[mc_idx] = self.fp_max_of_estar_interp(self.estar[mc_idx])
                self.fp_of_e_interp[mc_idx] = lambda e, mc_idx=mc_idx: self.fp_of_e_for_estar_interp(
                    e, self.estar[mc_idx])
                self.e_of_fp_interp[mc_idx] = lambda fp, mc_idx=mc_idx: self.e_of_fp_for_estar_interp(
                    fp, self.estar[mc_idx])
            print("Finished...................................\n") 


    def solve_snr(
        self,
        t=10. 
    ):
        print("Evaluating SNR...")  
        self.r_of_SNR8 = [0]*(len(self.mc_values))

        for mc_idx, mc_val in enumerate(self.mc_values): 
            self.r_of_SNR8[mc_idx] = lambda fp, mc_idx=mc_idx: ed.roffmSNR8(
                100., 
                fp, 
                self.e_of_fp_interp[mc_idx](fp),
                self.m1_values[mc_idx]*ed.msun_in_kg, 
                self.m2_values[mc_idx]*ed.msun_in_kg, 
                t
            )/((10.**6)*ed.pc_in_meters)

    def initialize_fp(self, fpmin_forced=1.0e-9):
        if (self.e_of_fp_interp==None):
            print("Need to 'solve_evolution(...)' first...")
            return

        print("Initializing binary peak frequencies...") 
        self.fp0 = [0]*(len(self.mc_values))
    
        for mc_idx, mc_val in enumerate(self.mc_values):
            self.fp0[mc_idx] = ed.dtdfp_sampler(
                mc_val, 
                self.e_of_fp_interp[mc_idx], 
                max(self.fp_min[mc_idx], fpmin_forced), 
                #self.fp_max[mc_idx],
                10.,  
                'log'
            ) 

    def solve_theta(self, e_cut, verbose=0):
        print("Evaluating SNR and eccentricity template cutoffs...") 
        self.e_cut = e_cut
        self.theta_cut = [0]*(len(self.mc_values))
        self.fp_cut = [0]*(len(self.mc_values))
        self.t_merge = [0]*(len(self.mc_values))
        self.ligo_cut = [0]*(len(self.mc_values))

        for mc_idx, mc_val in enumerate(self.mc_values): 
            self.fp_cut[mc_idx] = max(self.fp_of_e_interp[mc_idx](e_cut), self.fp_min[mc_idx])

            rsnr8 = self.r_of_SNR8[mc_idx](self.fp0[mc_idx])

            if (self.fp0[mc_idx]<self.fp_cut[mc_idx]): 
                self.theta_cut[mc_idx]=0
            elif (rsnr8 < self.chi_values[mc_idx]):
                self.theta_cut[mc_idx]=0
            else: 
                self.theta_cut[mc_idx]=1
            
            self.t_merge[mc_idx] = ed.tmerge(
                self.fp0[mc_idx], 
                self.e_of_fp_interp[mc_idx](self.fp0[mc_idx]), 
                self.m1_values[mc_idx]*ed.msun_in_kg, 
                self.m2_values[mc_idx]*ed.msun_in_kg
            ) 

            if (self.t_merge[mc_idx]<(10.*ed.year_in_seconds)): 
                self.ligo_cut[mc_idx]=1  

            if (verbose>0): 
                print(
                    "Binary: "+f"{mc_idx:d}"+", " 
                    +"fp0: "+f"{np.log10(self.fp0[mc_idx]):.3e}"+", "
                    +"log10(fpcut): "+f"{np.log10(self.fp_cut[mc_idx]):.3e}"+", "
                    +"r(snr>8): "+f"{rsnr8:.2f}"+", "
                    +"chi: "+f"{self.chi_values[mc_idx]:.2f}"+", "
                    +"theta_cut: "+f"{self.theta_cut[mc_idx]:d}"
                    +"t_merge: "+f"{self.t_merge[mc_idx]/ed.year_in_seconds:.2f}"+" yrs"
                    +"ligo_cut: "+f"{self.ligo_cut[mc_idx]:d}"
                )


    #def count_N(self, log10fmin, log10fmax, Max_N_binaries=-1):
    #    totalcount = 0.
    #    for mc_idx, mc_val in enumerate(self.mc_values[0:Max_N_binaries]):
    #        integrand = lambda log10fp: (
    #            1.
    #            * np.power(np.power(10., log10fp), -11./3.)
    #            / np.power(np.power(10., -1.5), -11./3.)
    #            * ed.F_script(self.e_of_fp_interp[mc_idx](np.power(10., log10fp)))
    #            /4.383
    #        )
    #        integral = scipy.integrate.quad(integrand, log10fmin, log10fmax)[0]
    #        counts = integral * self.theta_cut[mc_idx](np.power(10., log10fmin))
    #        totalcount += counts
    #        #print(integral, ", ", counts, ", ", totalcount)
    #    return totalcount

    #def count_N(self, log10fmin, log10fmax, binary_subset_count=None): 
    #    self.N_counts = np.zeros(len(self.mc_values))

    #    if (binary_subset_count==None): 
    #        nmax = len(self.mc_values)-1
    #    else: 
    #        nmax = int(binary_subset_count) 

    #    for mc_idx, mc_val in enumerate(self.mc_values[0:nmax]):
    #        if (log10fmax < np.log10(self.fp_cut[mc_idx])):
    #            self.N_counts[mc_idx]=0.
    #        else: 
    #            integrand = lambda log10fp, mc_idx=mc_idx: (
    #                ed.dtdfp(
    #                    self.mc_values[mc_idx], 
    #                    np.power(10., log10fp), 
    #                    self.e_of_fp_interp[mc_idx](np.power(10., log10fp))
    #                ) 
    #                * self.theta_cut[mc_idx](np.power(10., log10fp))
    #            )


    #            self.N_counts[mc_idx] = scipy.integrate.quad(integrand, log10fmin, log10fmax)[0]
    #            print("... ", self.N_counts[mc_idx], " for log10fpmin = ", log10fmin)            

    def count_N(self, log10fmin, log10fmax): 
        self.N_counts = 0

        for mc_idx, mc_val in enumerate(self.mc_values): 
            log10fp = np.log10(self.fp[mc_idx])

            if ((log10fp > log10fmin) and (log10fp < log10fmax)):
                self.N_counts += self.theta_cut[mc_idx]

        return self.N_counts 

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
            + 'Data columns: chi, m1, m2, mc, estar'
        )
        
        data = np.stack(
            (
                self.chi_values, 
                self.m1_values, 
                self.m2_values, 
                self.mc_values,
                self.estar
            ),
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
        w.estar = data[:, 4]
        w.world_mass = np.sum(w.m1_values + w.m2_values)

        w.redshifts = np.zeros(len(w.estar))
        cosmo = ed.Cosmology()
        for idx, val in enumerate(w.m1_values): 
            chi = w.chi_values[idx]
            w.redshifts[idx] = cosmo.redshift_solver(chi, 0.1) 
            w.chi_values[idx] = chi/(1.+w.redshifts[idx])
            w.mc_values[idx] = w.redshifts[idx]*w.mc_values[idx] 

        print("World loaded: ")
        print("\t Max comoving distance [Mpc]: ", chi_max)
        print("\t Annual Merger Rate [yr^-1 Mpc^-3]: ", annual_merger_rate)
        print("\t Mass distribution: ", mass_distribution)
        print("\t Mass args: ", mass_args)

        return w
      
    else: 
        print("Invalid filepath...")
        

