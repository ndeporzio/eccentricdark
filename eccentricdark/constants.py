import eccentricdark as ed

G = 6.674e-11 #Units: N * m^2 * kg^-2
c = 3.0e8 #Units: m * s^-1 
year_in_seconds = 31536000.
au_in_meters = 1.496e11
pc_in_meters = 3.086e16
msun_in_kg = 2.0e30

H0 = 67.7 # Units km/s/Mpc
omega_r = (10.**-4.)
omega_m = 0.311*(0.674**2.)
omega_k = 0.
omega_lambda = 0.689*(0.674**2.)
h = H0/100.
Omega_r = omega_r/(h**2.)
Omega_m = omega_m/(h**2.)
Omega_k = omega_k/(h**2.)
Omega_lambda = omega_lambda/(h**2.)

gamma = 1.1954 # 1907.02283

dtevolve = ((1.0)*year_in_seconds)
dt_rescale = 0.99 #factor to multiply dt by when a[t] or fp[t] begin to change to quickly
time_integral_cutoff_factor = 0.2 #when dt step equals this factor times initial dt, return value of SNR up to this point. Purpose is to cutoff the near merger region when dt gets very small and integrand gets very large and we're adding the produce of small and large numbers. 
invcdf_density_factor = 1000

e_interp_offset_default = 1.0e-9
# Can't interpolate e[fp] exactly at e=0 and e=1
# so we have to offset from these values by some
# amount in order to build an interpolation table. 
# This parameter sets what that offset should be. 

e_bin_count_default = int(1e6)
# When building e[fp] interpolation tables, this
# is how many values of e to calculate between 
# 0 and 1 to build interpolation table

lookup_bbh_evolution = False
