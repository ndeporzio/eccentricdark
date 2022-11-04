# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 14:22:26 2022

@author: Red October
"""

import eccentricdark as ed
import matplotlib.pyplot as plt
import numpy as np
import scipy
import os

fpbins = np.logspace(-2.7, -1.6, 11)


#isolated BBH

print('\n\n------------------------------------')
print('Begining Isolated BBH Calculation\n')
print('------------------------------------\n\n')

world_isolated = ed.World(chi_max=100., annual_merger_rate=50.0e-7)

world_isolated.populate_world(
    mass_distribution='1907.02283', 
    chi_distribution='uniform', 
    chi_args=[world_isolated.chi_max],
    cosmo=False)


world_isolated.initialize_eccentricities('isolated')
world_isolated.solve_evolution(mode="interpolate")

#initialize_fp checked to be good
world_isolated.initialize_fp(fpmin_forced=10**(-2.7))
world_isolated.solve_snr()
world_isolated.solve_theta(e_cut=0.9, verbose=1)
counts_isolated = np.zeros((11,))


for i in range(len(fpbins)-1):
    counts_isolated[i] = world_isolated.count_N(fpbins[i], fpbins[i+1])

world_isolated.save('./Isolated_World.txt')
    
plt.hist(np.log10(fpbins), weights=counts_isolated, density=True, histtype='step', color='blue', label='Isolated')

#ejected BBH

#print('\n\n------------------------------------')
#print('Begining Ejected BBH Calculation\n')
#print('------------------------------------\n\n')


#world_ejected = ed.World(chi_max=100., annual_merger_rate=50.0e-5)

#world_ejected.populate_world(
#    mass_distribution='1907.02283', 
#    chi_distribution='uniform', 
#     chi_args=[world_ejected.chi_max],
#     cosmo=False)

# world_ejected.initialize_eccentricities('ejected')
# world_ejected.solve_evolution(mode="interpolate")
# world_ejected.initialize_fp(fpmin_forced=10**(-2.7))
# world_ejected.solve_snr()
# world_ejected.solve_theta(e_cut=0.9, verbose=1)
# counts_ejected = np.zeros((11,))


# for i in range(len(fpbins)-1):
#     counts_ejected[i] = world_ejected.count_N(fpbins[i], fpbins[i+1])
   
# world_ejected.save('./Ejected_World.txt')
    
# plt.hist(np.log10(fpbins), weights=counts_ejected, density=True, histtype='step', color='red', label='Ejected')

# #in-cluster BBH

# print('\n\n------------------------------------')
# print('Begining In-Cluster BBH Calculation\n')
# print('------------------------------------\n\n')

# world_incluster = ed.World(chi_max=100., annual_merger_rate=50.0e-5)

# world_incluster.populate_world(
#     mass_distribution='1907.02283', 
#     chi_distribution='uniform', 
#     chi_args=[world_incluster.chi_max],
#     cosmo=False)

# world_incluster.initialize_eccentricities('incluster')
# world_incluster.solve_evolution(mode="interpolate")
# world_incluster.initialize_fp(fpmin_forced=10**(-2.7))
# world_incluster.solve_snr()
# world_incluster.solve_theta(e_cut=0.9, verbose=1)
# counts_incluster = np.zeros((11,))


# for i in range(len(fpbins)-1):
#     counts_incluster[i] = world_incluster.count_N(fpbins[i], fpbins[i+1])
    
# world_incluster.save('./World_Incluster.txt')    
    
# plt.hist(np.log10(fpbins), weights=counts_incluster, density=True, histtype='step', color='orange', label='In Cluster')

# #galcenter

# print('\n\n------------------------------------')
# print('Begining Gal Center BBH Calculation\n')
# print('------------------------------------\n\n')

# world_galcenter = ed.World(chi_max=100., annual_merger_rate=50.0e-5)

# world_galcenter.populate_world(
#     mass_distribution='1907.02283', 
#     chi_distribution='uniform', 
#     chi_args=[world_galcenter.chi_max],
#     cosmo=False)

# world_galcenter.initialize_eccentricities('galcenter')
# world_galcenter.solve_evolution(mode="interpolate")
# world_galcenter.initialize_fp(fpmin_forced=10**(-2.7))
# world_galcenter.solve_snr()
# world_galcenter.solve_theta(e_cut=0.9, verbose=1)
# counts_galcenter = np.zeros((11,))


# for i in range(len(fpbins)-1):
#     counts_galcenter[i] = world_galcenter.count_N(fpbins[i], fpbins[i+1])
    
# world_galcenter.save('./GalCenter_World.txt')    
    
# plt.hist(np.log10(fpbins), weights=counts_galcenter, density=True, histtype='step', color='green', label='Gal Center')

# plt.legend()
# savepath = '.'
# plt.savefig(os.path.join(savepath, 'Figure_1.png'))
# plt.show()