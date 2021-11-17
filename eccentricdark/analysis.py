import scipy 
import numpy as np
import eccentricdark as ed


class Analysis: 

    def __init__(
        self,
        world, 
        estar, 
        fpstar,
        log10fpmin, 
        log10fpmax, 
        fpsteps, 
    ): 

        self.world = world
        self.estar = estar 
        self.fpstar = fpstar

        self.fpsteps = fpsteps

        self.fpmin = None
        self.fpmax = None


        fpcheck = np.logspace(log10fpmin, log10fpmax, 1001)

        print("Identifying lower fp bound...") 
        counter = 0
        while self.fpmin==None: 
            if (ed.e_solver(fpcheck[counter], self.fpstar, self.estar) != 1):
                self.fpmin = fpcheck[counter]
            else: 
                counter+=1

        print("Identifying upper fp bound...") 
        counter = 1000
        while self.fpmax==None: 
            if (ed.e_solver(fpcheck[counter], self.fpstar, self.estar) != 0):
                self.fpmax = fpcheck[counter]
            else:
                counter-=1 

        print("Initializing frequencies...") 
        self.fptable = np.geomspace(self.fpmin, self.fpmax, self.fpsteps) 

        self.dNdfp = np.zeros((len(self.world.mc_values), len(self.fptable)))
        self.N = np.zeros((len(self.world.mc_values), len(self.fptable)))
        for mcidx, mcval in enumerate(self.world.mc_values):
            if ((mcidx%int(np.floor(len(self.world.mc_values)/100.)))==0):
                print(
                    'Sample evaluation is '
                    + str(np.floor(100.*mcidx/len(self.world.mc_values)))
                    + '% complete...'
                )
            for fpidx, fpval in enumerate(self.fptable):
                self.N[mcidx, fpidx] = ed.N_estarfixed(
                    R=self.world.annual_merger_rate, 
                    f=fpval,
                    chimax=self.world.chi_max, 
                    mass_distribution=self.world.mass_distribution, 
                    estar=self.estar, 
                    fpstar=self.fpstar
                )
                #if mcidx==0:  
                #    self.dNdfp[mcidx, fpidx] = ed.dNdfp_estarfixed_mcfixed(
                #        R=self.world.annual_merger_rate, 
                #        chimax=self.world.chi_max, 
                #        mc=mcval, 
                #        fp=fpval, 
                #        fpstar=self.fpstar, 
                #        estar=self.estar)
                #else: 
                #    self.dNdfp[mcidx, fpidx] = (
                #        self.dNdfp[0, fpidx] 
                #        * np.power(mcval/self.world.mc_values[0], -5./3.)
                #    )


        self.dNdfp_interp = [
            scipy.interpolate.interp1d(self.fptable, self.dNdfp[idx, :])
            for idx, val in enumerate(self.fptable)]

        self.dNdfp_normalization = scipy.integrate.quad(
            self.dNdfp_interp[0], 
            self.fptable[0], 
            self.fptable[-1], 
            limit=1000
        )[0]        

        self.dNdfp_normed = self.dNdfp/self.dNdfp_normalization

        self.dNdfp_normed2 = np.array([
            self.dNdfp[midx][0:-1]/np.sum(self.dNdfp[midx][0:-1]*np.diff(self.fptable))            
        for midx, mval in enumerate(self.world.mc_values)])

        self.dNdfp_normed_interp = [
            scipy.interpolate.interp1d(
                self.fptable[0:-1], 
                self.dNdfp_normed2[idx, :],
                bounds_error=False,
                fill_value="extrapolate")
            for idx, val in enumerate(self.fptable)]

    def fpi_cdf(self, mass, fp): 
        massidx = np.argmin(np.abs(self.world.mc_values-mass))
        integrand = lambda x: self.dNdfp_normed_interp[massidx](x)
        return scipy.integrate.quad(integrand, self.fptable[0], fp, limit=1000)[0]

    def sample_fpi(self, mass):
        rand = np.random.random() 
        searchfunc = lambda f: self.fpi_cdf(mass, f) - rand
        try: 
            fpi = scipy.optimize.bisect(
                searchfunc, 
                (1.-1e-2)*self.fptable[0], 
                (1.+1e-2)*self.fptable[-1]
            ) 
            return fpi
        except: 
            return 0

    def generate_fi(self): 
        self.fpi = np.zeros(len(self.world.mc_values))
        for midx, mval in enumerate(self.world.mc_values): 
            self.fpi[midx] = self.sample_fpi(mval)
        if ((midx%int(np.floor(len(self.world.mc_values)/100.)))==0):
           print(
               'Frequency initialization is '
               + str(np.floor(100.*midx/len(self.world.mc_values)))
               + '% complete... fpi = ' + f'{self.fpi[midx]:.3e}'
           )

    def sample_fpi_2(self, mass, massidx=None):
        if (massidx==None): 
            massidx = np.argmin(np.abs(self.world.mc_values-mass)) 
        rand = np.random.random()
        cdf = self.dNdfp_normed2[massidx, :]*np.diff(self.fptable)
        for fidx, fval in enumerate(self.fptable): 
            if np.sum(cdf[0:fidx])>rand: 
                return fval 
    
    def generate_fi_2(self): 
        self.fpi = np.zeros(len(self.world.mc_values))
        for midx, mval in enumerate(self.world.mc_values):
            self.fpi[midx] = self.sample_fpi_2(mval)
            if ((midx%int(np.floor(len(self.world.mc_values)/100.)))==0):
                print(
                    'Frequency initialization is '
                    + str(np.floor(100.*midx/len(self.world.mc_values)))
                    + '% complete... fpi = ' + f'{self.fpi[midx]:.3e}'
                )
    







