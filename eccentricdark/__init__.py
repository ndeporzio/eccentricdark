# The lines in this file are executed upon importing the 
# 'eccentricdark' package in a Python environment.

from .constants import (
    G, 
    c, 
    year_in_seconds, 
    au_in_meters, 
    pc_in_meters, 
    msun_in_kg, 
    H0, 
    omega_r, 
    omega_m, 
    omega_k, 
    omega_lambda, 
    h, 
    Omega_r, 
    Omega_m, 
    Omega_k, 
    Omega_lambda,
    gamma
)
from .world import (
    World, 
    load_world
)
from .massdistributions import (
    mass_distribution_sampler
)
from .cosmology import (
    Cosmology,
    load_cosmology,
    H, 
    comoving_distance, 
    comoving_volume, 
    redshift_solver
)
from .equations import (
    m_chirp, 
    G_script, 
    H_script, 
    e_solver, 
    F_script, 
    dtdfp, 
    dNdfp_estarfixed_mcfixed
)
from .analysis import (
    Analysis
)

#from ./lookup/
#from .snr import ()
#from .plotting import ()

print(
    "You have loaded the eccentricdark - a Python package for modeling"
    + " experimental sensitivity to black hole binary formation "
    + " channels through eccentricity observations... \n"
    + "See 'Tutorial.ipynb' for use instructions..." 
)
