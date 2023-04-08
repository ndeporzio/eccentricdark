import copy 
import eccentricdark as ed
import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import sys
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib import rc

#sns.set()
rc('font', **{'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
matplotlib.rcParams.update({
    "font.weight" : "bold",
    "font.size" : 30,
    "axes.labelsize" : 60,
    "axes.labelpad" : 10.0,
    "xtick.labelsize" : 30,
    "xtick.major.size" : 30,
    "xtick.major.width" : 3,
    "xtick.major.pad" : 15,
    "xtick.minor.size" : 20,
    "xtick.minor.width" : 2,
    "xtick.direction" : "in",
    "ytick.labelsize" : 30,
    "ytick.major.size" : 30,
    "ytick.major.width" : 3,
    "ytick.major.pad" : 10,
    "ytick.minor.size" : 20,
    "ytick.minor.width" : 2,
    "ytick.direction" : "in",
    "legend.fontsize" : 60,
    "figure.dpi" : 100,
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True,
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

imbh = False
fig1 = False 
fig2 = True 
fig3 = False 
fig4 = False
extra1 = False
extra2 = False 

#Collect output directory 
if(len(sys.argv)!=2): 
    sys.exit("\nERROR: Please (only) specify the output directory when calling this script!\n") 
else: 
    savepath = sys.argv[1]
    print("\nSaving results to: ", savepath, "\n")

#I think this is from IMBH plots? 
if imbh==True: 
    fp = np.logspace(-4, 1, 51)
    m = np.geomspace(1.0e1*ed.msun_in_kg, 1.0e6*ed.msun_in_kg, 6)
    y = np.array([[
            ed.roffmSNR8(100., f, 0., mval, mval/4., 10.)
            for f in fp]
            for mval in m]
        )
    y2 = np.nan_to_num(y)
    roffmSNR8intpl = scipy.interpolate.interp2d(fp, m, y2, kind='linear')
    sns.set_palette("magma", len(m))
    fp = np.logspace(-4, 1, 601)
    plt.figure(figsize=(15, 7.5))
    for mass in m:
        plt.plot(
            np.log10(fp),
            np.log10(roffmSNR8intpl(fp, mass)/((10.**6)*ed.pc_in_meters)),
            label="M = "+f'{mass/ed.msun_in_kg:.2e}'+r' $M_{sun}$'
        )

    ax = plt.gca()
    lisa_rect = Rectangle(
        (-2.6, -3.),
        1.0,
        9.,
        facecolor='blue',
        alpha=0.2)
    tian_rect = Rectangle(
        (-1.8, -3.),
        1.0,
        9.,
        facecolor='green',
        alpha=0.2
    )
    decigo_rect = Rectangle(
        (-1., -3.),
        2.,
        9.,
        facecolor='red',
        alpha=0.2
    )
    ligo_rect = Rectangle(
        (0.477, -3),
        1,
        9.,
        facecolor='grey',
        alpha=0.2
    )
    ax.add_patch(lisa_rect)
    ax.add_patch(tian_rect)
    ax.add_patch(decigo_rect)
    ax.add_patch(ligo_rect)

    plt.text(-2.5, -2.,"LISA")
    plt.text(-1.5, -2.,"TianQin")
    plt.text(-0.5, -2.,"DECIGO")
    plt.text(0.6, -2.,"LIGO")
    plt.ylim((-3., 6.))
    plt.xlim((-4, 1.))
    plt.grid(True, which='both', axis='both')
    plt.xlabel(r"log($f_p/$[Hz])", fontsize=20)
    plt.ylabel(r"log(r/[Mpc]) s.t. SNR>8", fontsize=20)
    plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
    plt.legend()

    plt.savefig(os.path.join(savepath, 'Figure_1.png'))

#Check LISA SNR - Figure 1 
if fig1==True:
    print("Generating 1907.02283 Figure 1...") 
    fp_table = np.logspace(-4, 2, 300)
    estar_table = np.array([1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6])

    e_of_fp = np.zeros((len(estar_table), len(fp_table)))
    for estaridx, estarval in enumerate(estar_table):
        print("\tEvaluating e(fp) for estar = ", estarval)
        e_of_fp_interp = ed.e_to_fp_interpolator(
            fp_star=ed.fp_star_default,
            e_star=estarval
        )[3] 
        for fpidx, fpval in enumerate(fp_table): 
            e_of_fp[estaridx, fpidx] = e_of_fp_interp(fpval)

    colors=sns.dark_palette("#69d", len(estar_table), reverse=False)

    plt.figure(figsize=(15, 7.5))
    plt.title("Fig. 1 in 1907.02283", fontsize=20)
    plt.xlabel(r"$f_p$/Hz", fontsize=20)
    plt.ylabel(r"$e$", fontsize=20)
    for estaridx, estarval in enumerate(estar_table): 
        plt.plot(fp_table, e_of_fp[estaridx, :], label=r"$e_*=$"+f'{estarval:.1e}', color=colors[estaridx])
    plt.plot([1.0e-4, 1.0e2], [1., 1.], color='black', label=r"e=1.0")
    plt.plot([1.0e-4, 1.0e2], [2.0e-3, 2.0e-3], color='black', alpha=0.2)
    plt.xlim((1.0e-4, 1.0e2))
    plt.ylim((1.0e-3, 2.))
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(axis='both', which='both')

    ymin=1.0e-3
    ymax=2.0
    lisa_fmin=np.power(10., -2.6)
    lisa_fmax=np.power(10., -1.6)
    tian_fmin=np.power(10., -1.8)
    tian_fmax=np.power(10., -0.8)
    decigo_fmin=np.power(10., -1.0)
    decigo_fmax=np.power(10., 1.0)
    ligo_fmin=3.
    ligo_fmax=200. 

    ax = plt.gca()
    lisa_rect = Rectangle(
        (lisa_fmin, ymin), #Bottom-left corner
        lisa_fmax-lisa_fmin, #Width
        ymax-ymin, #Height
        facecolor='grey',
        alpha=0.2)
    tian_rect = Rectangle(
        (tian_fmin, ymin),
        tian_fmax-tian_fmin,
        ymax-ymin,
        facecolor='green',
        alpha=0.2
    )
    decigo_rect = Rectangle(
        (decigo_fmin, ymin),
        decigo_fmax-decigo_fmin,
        ymax-ymin,
        facecolor='blue',
        alpha=0.2
    )
    ligo_rect = Rectangle(
        (ligo_fmin, ymin),
        ligo_fmax-ligo_fmin,
        ymax-ymin,
        facecolor='grey',
        alpha=0.2
    )
    ax.add_patch(lisa_rect)
    ax.add_patch(tian_rect)
    ax.add_patch(decigo_rect)
    ax.add_patch(ligo_rect)

    plt.text(lisa_fmin, 3.0e-2, "LISA", rotation='vertical', fontsize=20)
    plt.text(tian_fmin, 3.0e-2, "TianQin", rotation='vertical', fontsize=20)
    plt.text(decigo_fmin, 3.0e-2, "DECIGO", rotation='vertical', fontsize=20)
    plt.text(ligo_fmin, 3.0e-2, "LIGO", rotation='vertical', fontsize=20)

    estar_merge = np.logspace(-8, -1, 50)
    t1month = np.zeros((len(estar_merge), 2))
    t12month = np.zeros((len(estar_merge), 2))
    t120month = np.zeros((len(estar_merge), 2))
    t1200month = np.zeros((len(estar_merge), 2))
    
    for e_idx, e_val in enumerate(estar_merge):
        fpmin, fpmax, fp_of_e, e_of_fp = ed.e_to_fp_interpolator(fp_star=ed.fp_star_default, e_star=e_val)
        
        print("\tBuilding merge time : ", e_idx+1, " of ", len(estar_merge))
        tmerge_1month = lambda fp: ed.tmerge(fp, e_of_fp(fp), 30.*ed.msun_in_kg, 30.*ed.msun_in_kg)-(1./12.)*ed.year_in_seconds
        tmerge_12month = lambda fp: ed.tmerge(fp, e_of_fp(fp), 30.*ed.msun_in_kg, 30.*ed.msun_in_kg)-(1.)*ed.year_in_seconds
        tmerge_120month = lambda fp: ed.tmerge(fp, e_of_fp(fp), 30.*ed.msun_in_kg, 30.*ed.msun_in_kg)-(10.)*ed.year_in_seconds
        tmerge_1200month = lambda fp: ed.tmerge(fp, e_of_fp(fp), 30.*ed.msun_in_kg, 30.*ed.msun_in_kg)-(100.)*ed.year_in_seconds
    
        interp_lims = ed.e_to_fp_interpolator(fp_star=ed.fp_star_default, e_star=e_val)
    
        if (t1month[e_idx-1, 1] < (1.-1.0e-3)):
            fp_1month = scipy.optimize.bisect(tmerge_1month, interp_lims[0], 1.0e0, xtol=1.0e-8, rtol=1.0e-8)
        if (t12month[e_idx-1, 1] < (1.-1.0e-3)):
            fp_12month = scipy.optimize.bisect(tmerge_12month, interp_lims[0], 1.0e0, xtol=1.0e-8, rtol=1.0e-8)   
        if (t120month[e_idx-1, 1] < (1.-1.0e-3)):
            fp_120month = scipy.optimize.bisect(tmerge_120month, interp_lims[0], 1.0e0, xtol=1.0e-8, rtol=1.0e-8)   
        if (t1200month[e_idx-1, 1] < (1.-1.0e-3)):
            fp_1200month = scipy.optimize.bisect(tmerge_1200month, interp_lims[0], 1.0e0, xtol=1.0e-8, rtol=1.0e-8)
            
        t1month[e_idx, 0] = fp_1month
        t12month[e_idx, 0] = fp_12month
        t120month[e_idx, 0] = fp_120month
        t1200month[e_idx, 0] = fp_1200month
    
        t1month[e_idx, 1] = e_of_fp(fp_1month)
        t12month[e_idx, 1] = e_of_fp(fp_12month)
        t120month[e_idx, 1] = e_of_fp(fp_120month)
        t1200month[e_idx, 1] = e_of_fp(fp_1200month)


    colors=sns.color_palette("Reds", 4)
    plt.plot(t1month[:,0], t1month[:,1], color=colors[0], label="1 month", linestyle='dashed') 
    plt.plot(t12month[:,0], t12month[:,1], color=colors[1], label="1 year", linestyle='dashed') 
    plt.plot(t120month[:,0], t120month[:,1], color=colors[2], label="10 year", linestyle='dashed') 
    plt.plot(t1200month[:,0], t1200month[:,1], color=colors[3], label="100 year", linestyle='dashed') 

    plt.legend(fontsize=20)
    plt.savefig(os.path.join(savepath, 'Figure_1.png'))


# Check primordial distribution functions - Figure 2
if fig2==True: 
    plt.figure(figsize=(26.25, 7.5))

    # Sample the True Distributions and Plot
    nsample = 100000
    edges = np.logspace(-7, -2, 51)

    isolated_sample = ed.estar_sampler("isolated")(np.random.random(nsample))
    ejected_sample = ed.estar_sampler("ejected")(np.random.random(nsample))
    incluster_sample = ed.estar_sampler("incluster")(np.random.random(nsample))
    galcenter_sample = ed.estar_sampler("galcenter")(np.random.random(nsample))

    isolated_counts, bins = np.histogram(isolated_sample, edges)
    ejected_counts, bins = np.histogram(ejected_sample, edges)
    incluster_counts, bins = np.histogram(incluster_sample, edges)
    galcenter_counts, bins = np.histogram(galcenter_sample, edges)

    plt.hist(bins[:-1], bins, weights=isolated_counts/max(isolated_counts), 
        label=r'isolated', linewidth=2, histtype='stepfilled', color='grey', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=ejected_counts/max(ejected_counts), 
        label=r'ejected', linewidth=2, histtype='stepfilled', color='red', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=incluster_counts/max(incluster_counts), 
        label=r'in-cluster', linewidth=2, histtype='stepfilled', color='purple', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=galcenter_counts/max(galcenter_counts), 
        label=r'gal. center', linewidth=2, histtype='stepfilled', color='blue', alpha=0.5)

    #Sample a gaussian and plot
    def gauss(x, u=1.0e-3, sig=3.0e-4): 
        return (1./(sig*np.sqrt(2.*np.pi)))*np.power(np.e, (-1./2.)*np.power((x-u)/sig, 2.))

    gauss_sample = ed.generate_invcdf(gauss, 1.0e-4, 1.0e-2, 'log')(np.random.random(nsample))
    gauss_counts, bins = np.histogram(gauss_sample, edges)
    plt.hist(bins[:-1], bins, weights=gauss_counts/max(gauss_counts),
        label=r'Gaussian', linewidth=2, histtype='stepfilled', color='green', alpha=0.5)    
    plt.plot(np.logspace(-5, -2, 301), gauss(np.logspace(-5, -2, 301))/max(gauss(np.logspace(-5, -2, 301))), 
        linestyle='dashed', color='green', linewidth=4)

    # Plotting True Distributions
    plt.plot(np.power(10., ed.fieldData[:,0]), ed.fieldData[:,1]/max(ed.fieldData[:,1]), 
        linewidth=4, color="grey")
    plt.plot(np.power(10., ed.ejectedData[:,0]), ed.ejectedData[:,1]/max(ed.ejectedData[:,1]), 
        linewidth=4, color="red")
    plt.plot(np.power(10., ed.inclusterData[:,0]), ed.inclusterData[:,1]/max(ed.inclusterData[:,1]), 
        linewidth=4, color="purple")
    plt.plot(np.power(10., ed.galcenterData[:,0]), ed.galcenterData[:,1]/max(ed.galcenterData[:,1]), 
        linewidth=4, color="blue")


    # Plotting Log-Normal Approximations to True Distributions
    xplot = np.logspace(-7, -4, 200)
    fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss,
            np.power(10., ed.fieldData[:,0]),
            ed.fieldData[:,1]/max(ed.fieldData[:,1]),
            [-12.4, -0.43, -0.10, -1.4e-2, 1.6e-6]
    )
    plt.plot(xplot, ed.multigauss(
            y = xplot,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        ), color="grey", linestyle='dashed', linewidth=4)

    xplot = np.logspace(-7, -3, 200)
    fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss,
            np.power(10., ed.ejectedData[:,0]),
            ed.ejectedData[:,1]/max(ed.ejectedData[:,1]),
            [-8.8, -2.1e-1, -5.0e-2, -5.2e-3, 9.1e-6]
    )
    plt.plot(xplot, ed.multigauss(
            y = xplot,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        ), color="red", linestyle='dashed', linewidth=4)

    xplot = np.logspace(-7, -2, 200)
    fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss,
            np.power(10., ed.inclusterData[:,0]),
            ed.inclusterData[:,1]/max(ed.inclusterData[:,1]),
            [-7.2, -3.9e-1, -1.2e-1, -1.7e-2, 1.6e-4]
    )
    plt.plot(xplot, ed.multigauss(
            y = xplot,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        ), color="purple", linestyle='dashed', linewidth=4)

    xplot = np.logspace(-7, -2, 200)
    fit1, fit2 = scipy.optimize.curve_fit(
            ed.multigauss,
            np.power(10., ed.galcenterData[:,0]),
            ed.galcenterData[:,1]/max(ed.galcenterData[:,1]),
            [-4.7, -5.5e-1, -1.8e-1, -1.8e-2, 1.0e-3]
    )
    plt.plot(xplot, ed.multigauss(
            y = xplot,
            u = fit1[0],
            a = fit1[1],
            b = fit1[2],
            c = fit1[3],
            k = fit1[4]
        ), color="blue", linestyle='dashed', linewidth=4)

    # Plotting
    plt.legend(fontsize=20)
    plt.ylabel(r"$e_0*$-distribution", fontsize=20)
    plt.xlabel(r"$e_0$ at $f_{p0}=10$ Hz (i.e. $e_*$ for $f_{p*}=10$ Hz", fontsize=20)
    plt.xscale('log') 
    plt.xlim((1.0e-7, 1.0e-2))
    plt.savefig(os.path.join(savepath, 'Figure_2.png'))


# Check LISA SNR - Figure 3
if fig3==True:
    #Fig 3 left 
    xplot = np.linspace(0., 1., 501)
    snr_ref = ed.SNR_LISA(
        r = 1.0e8,
        fp0 = 5.0e-3, 
        e0 = 1.0e-12,  
        m1 = 50., 
        m2 = 50., 
        ttoday = 10. * ed.year_in_seconds
    )
    snr_vals = np.array([ed.SNR_LISA(
        r = 1.0e8,
        fp0 = 5.0e-3, 
        e0 = xval,  
        m1 = 50., 
        m2 = 50., 
        ttoday = 10. * ed.year_in_seconds
    ) for xidx, xval in enumerate(xplot)])
    yplot = snr_vals/snr_ref
    
    plt.figure(figsize=(7.5,7.5))
    plt.plot(xplot, np.nan_to_num(yplot))
    plt.xlabel("e", fontsize=20)
    plt.ylabel("SNR(e)/SNR(e=0)", fontsize=20)
    plt.grid(True, which='both', axis='both')
    plt.xlim((0.,1.02))
    plt.ylim((0.,1.05))
    plt.savefig(os.path.join(savepath, 'Figure_3a.png'))
    
    # Fig 3 right
    xplot = np.linspace(0., 1., 501)
    
    def Nfactor(e): 
        return (np.power(1.-e, 9./4.)*ed.F_script(e))
    
    N_ref = Nfactor(e=0.)
    N_vals = Nfactor(xplot)
    yplot = N_vals/N_ref
    
    plt.figure(figsize=(7.5,7.5))
    plt.plot(xplot, np.nan_to_num(yplot))
    plt.xlabel("e", fontsize=20)
    plt.ylabel("N(e)/N(e=0)", fontsize=20)
    plt.grid(True, which='both', axis='both')
    plt.xlim((0.,1.02))
    plt.ylim((0.,1.65))
    plt.savefig(os.path.join(savepath, 'Figure_3b.png'))


# Figure 4
if fig4==True: 
    #Setup world size
    max_chi_in_Mpc = 3000.
    world1 = ed.World(
        chi_max=max_chi_in_Mpc, #Max comoving radius of spherical world
        annual_merger_rate=50.0e-9 #Units: binaries/year/Mpc^3
    )
    
    #Populate world with binaries 
    world1.populate_world(
        mass_distribution='1907.02283', #Use same mass distribution 
        chi_distribution='uniform', #Fix binaries at one comoving distance 
        chi_args=[max_chi_in_Mpc], #Fix binaries at 100 Mpc comoving distance
        cosmo=False
    )
    
    #Visualize binary mass distribution
    plt.figure(figsize=(15, 7.5))
    plt.hist(world1.m1_values, density=True, bins=200, label=r'$m_1$ Samples')
    plt.hist(world1.m2_values, density=True, bins=200, label=r'$m_2$ Samples', histtype='step')
    plt.plot(
        np.arange(5, 50, 0.1), 
        (1.3 / (np.power(5., -1.3) - np.power(50., -1.3)))*np.power(np.arange(5, 50, 0.1), -2.3), 
        label=r'PDF $\propto m_1^{-2.3}$')
    plt.legend(fontsize=20)
    plt.xlabel('m1 Mass [Solar Mass]', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.grid(True, which='both', axis='both')
    plt.savefig(os.path.join(savepath, 'mass_dist_1.png'))
    
    plt.figure(figsize=(15, 7.5))
    plt.hist(world1.mc_values, density=True, bins=200, label=r'$m_c$ Samples')
    plt.hist(world1.m1_values, density=True, bins=200, label=r'$m_1$ Samples', histtype='step')
    plt.plot(np.arange(5, 50, 0.1), (1.3 / (np.power(5., -1.3) - np.power(50., -1.3)))*np.power(np.arange(5, 50, 0.1), -2.3), label=r'$m_1$ PDF')
    plt.legend(fontsize=20)
    plt.xlabel('mc Mass [Solar Mass]', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.grid(True, which='both', axis='both')
    plt.savefig(os.path.join(savepath, 'mass_dist_2.png'))
    
    #Visualize binary position distribution
    nbins = 200
    xplot = np.arange(0., world1.chi_max, world1.chi_max/nbins)
    expected = (4./3.)*np.pi*np.power(world1.chi_max, 3.)*world1.annual_merger_rate
    plt.figure(figsize=(15, 7.5))
    plt.hist(world1.chi_values, density=False, bins=nbins)
    plt.plot(
        xplot[0:-1], 
        (4.*np.pi*world1.annual_merger_rate)*np.power(xplot[0:-1], 2.)*np.diff(xplot),
        label=r'$\frac{dN}{d\chi} = 4 \pi \mathcal{R} \chi^2$'
    )
    plt.text(100, 20, "Total Binaries: "+f'{len(world1.chi_values):d}')
    plt.text(100, 15, "Expected Binaries: "+f'{expected:.1f}') 
    plt.xlabel(r'$\chi$ comoving [Mpc]', fontsize=20)
    plt.ylabel('Binary Counts', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid(True, which='both', axis='both')
    plt.savefig(os.path.join(savepath, 'chi_dist.png'))
    
    #Initialize binary eccentricities
    #world2 = copy.deepcopy(world1)
    #world3 = copy.deepcopy(world1)
    #world4 = copy.deepcopy(world1)
    #world5 = copy.deepcopy(world1)
    
    world1.initialize_eccentricities('fixed', np.power(10., -5.5))
    #world2.initialize_eccentricities('fixed', np.power(10., -4.0))
    #world3.initialize_eccentricities('fixed', np.power(10., -4.5))
    #world4.initialize_eccentricities('fixed', np.power(10., -5.0))
    #world5.initialize_eccentricities('fixed', np.power(10., -5.5))
    
    #Visualize binary eccentricity distributions
    plt.figure(figsize=(15, 7.5))
    plt.hist(world1.estar, density=False, bins=np.logspace(-6, -3, 7), label="world 1")
    #plt.hist(world2.estar, density=False, bins=np.logspace(-6, -3, 7), label="world 2")
    #plt.hist(world3.estar, density=False, bins=np.logspace(-6, -3, 7), label="world 3")
    #plt.hist(world4.estar, density=False, bins=np.logspace(-6, -3, 7), label="world 4")
    #plt.hist(world5.estar, density=False, bins=np.logspace(-6, -3, 7), label="world 5")
    plt.xlabel(r'$e_*$', fontsize=20)
    plt.ylabel('Counts', fontsize=20)
    plt.legend()
    plt.grid(True, which='both', axis='both')
    plt.xscale('log')
    plt.savefig(os.path.join(savepath, 'e_dist.png'))
    
    #Evaluate evolution e[fp|estar] for binaries
    world1.solve_evolution(mode="fixed")
    #world2.solve_evolution(mode="fixed")
    #world3.solve_evolution(mode="fixed")
    #world4.solve_evolution(mode="fixed")
    #world5.solve_evolution(mode="fixed")
    
    #Evaluate LISA SNR for each world
    print("Generating SNR functions...") 
    world1.solve_snr()
    #world2.solve_snr()
    #world3.solve_snr()
    #world4.solve_snr()
    #world5.solve_snr()
    
    #Count number of observable binaries given maximally observable eccentricity
    e_cutoffs=[0.01, 0.1, 0.4, 0.9]
    fpbins = np.logspace(-2.7, -1.5, 13)

    data1 = np.zeros((len(e_cutoffs), len(fpbins)-1))
    #data2 = np.zeros((len(e_cutoffs), len(fpbins)-1))
    #data3 = np.zeros((len(e_cutoffs), len(fpbins)-1))
    #data4 = np.zeros((len(e_cutoffs), len(fpbins)-1))
    #data5 = np.zeros((len(e_cutoffs), len(fpbins)-1))
    
    for ecutidx, ecutval in enumerate(e_cutoffs): 
        world1.solve_theta(e_cut=ecutval)
        #world2.solve_theta(e_cut=ecutval)
        #world3.solve_theta(e_cut=ecutval)
        #world4.solve_theta(e_cut=ecutval)
        #world5.solve_theta(e_cut=ecutval)
    
        for fp_idx, fp_val in enumerate(fpbins[0:-1]):
            world1.count_N2(
                np.log10(fpbins[fp_idx]), 
                np.log10(fpbins[fp_idx+1])
            )

            data1[ecutidx][fp_idx] = np.sum(world1.N_counts) 
        #    data2[ecutidx][fp_idx] = count_N(world2, np.log10(fp_val), np.log10(fp_val)+log10dfp)
        #    data3[ecutidx][fp_idx] = count_N(world3, np.log10(fp_val), np.log10(fp_val)+log10dfp)
        #    data4[ecutidx][fp_idx] = count_N(world4, np.log10(fp_val), np.log10(fp_val)+log10dfp)
        #    data5[ecutidx][fp_idx] = count_N(world5, np.log10(fp_val), np.log10(fp_val)+log10dfp)
    
    np.savetxt(savepath+"/count_N.txt", data1)


    np.savetxt(savepath+"/nocosmo_chi3000uniform_estar3p5.txt", data1)
    #np.savetxt(savepath+"/nocosmo_chi3000uniform_estar4p0.txt", data2)
    #np.savetxt(savepath+"/nocosmo_chi3000uniform_estar4p5.txt", data3)
    #np.savetxt(savepath+"/nocosmo_chi3000uniform_estar5p0.txt", data4)
    #np.savetxt(savepath+"/nocosmo_chi3000uniform_estar5p5.txt", data5)
    
    #Visualize observable number counts
    colors = [
        (193/255., 178/255., 240/255., 0.7), #purple  
        (178/255., 209/255., 240/255., 0.7), #blue 
        (178/255., 240/255., 209/255., 0.7), #green 
        (255/255., 224/255., 178/255., 0.7) #orange
    ]
    
    
#    plt.figure(figsize=(7.5, 7.5))
#    plt.hist(
#        np.log10(fpbins)[0:-1], 
#        weights=data5[0], 
#        bins=np.log10(fpbins), 
#        color=colors[0],  
#        label=r'$e < 0.01$', 
#        linewidth=4, 
#        histtype='step'
#    )
#    plt.hist(
#        np.log10(fpbins)[0:-1], 
#        weights=data5[1], 
#        bins=np.log10(fpbins), 
#        color=colors[1], 
#        label=r'$0.01 < e < 0.1$', 
#        linewidth=4, 
#        histtype='step'
#    )
#    plt.title("Fixed Formation $e_*$ Distribution, Static Cosmology", fontsize=20)
#    plt.legend(fontsize=20)
#    plt.xlabel(r"$\log (f_p/Hz)$", fontsize=20)
#    plt.ylabel(r"LISA N2A5 Observable Counts", fontsize=20)
#    plt.savefig(savepath+"/nocosmo_chivaried_estarfixed.png")


# Check LISA SNR Functions
if extra1==True: 
    x = np.logspace(-5, 0, 501)
    y = np.array([np.sqrt(ed.SnLISAdefault(xval)) for xval in x])
    plt.figure(figsize=(15, 7.5))
    plt.plot(x, y)
    plt.grid(True, which='both', axis='both')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$f_p$ [Hz]", fontsize=20)
    plt.ylabel(r"$(Noise)^{1/2}$", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
    plt.savefig(os.path.join(savepath, 'LISA_N2A5_1.png'))
    
    x = np.logspace(-5, 0, 501)
    y = np.array([np.power(xval, 0.5)*np.sqrt(ed.SnLISAdefault(xval)) for xval in x])
    plt.figure(figsize=(15, 7.5))
    plt.plot(x, y)
    plt.grid(True, which='both', axis='both')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$f_p$ [Hz]", fontsize=20)
    plt.ylabel(r"$(f_p \times Noise)^{1/2}$", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
    plt.savefig(os.path.join(savepath, 'LISA_N2A5_2.png'))
    
    fp = np.logspace(-4, 1, 51)
    m = np.arange(10.*ed.msun_in_kg, 100.*ed.msun_in_kg, 10.*ed.msun_in_kg)
    y = np.array([[
            ed.roffmSNR8(100., f, 0., mval, mval/4., 10.)
            for f in fp]
            for mval in m]
        )
    y2 = np.nan_to_num(y)
    roffmSNR8intpl = scipy.interpolate.interp2d(fp, m, y2, kind='linear')
    sns.set_palette("magma", len(m))
    fp = np.logspace(-4, 1, 501)
    plt.figure(figsize=(15, 7.5))
    for mass in m: 
        plt.plot(fp, roffmSNR8intpl(fp, mass)/((10.**6)*ed.pc_in_meters), label="M = "+f'{mass/ed.msun_in_kg:.2e}'+r' $M_{sun}$')
    plt.xscale('log')
    plt.grid(True, which='both', axis='both')
    plt.xlabel(r"$f_p$ [Hz]", fontsize=20)
    plt.ylabel(r"r [Mpc] s.t. SNR>8", fontsize=20)
    plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
    plt.legend()
    plt.savefig(os.path.join(savepath, 'LISA_N2A5_3.png'))

if extra2==True: 
    evolution = ed.e_to_fp_interpolator(10., np.power(10., -4.))
    e_of_fp_interp = evolution[3]
    fp_min = evolution[0]
    fp_max = evolution[1]
    mc = ed.m_chirp(10., 10.)

    xplot = np.logspace(np.log10(fp_min)+0.1, -1.5, 100)
    yplot = np.linspace(0., 1000., 1001)
    xx, yy = np.meshgrid(xplot, yplot)

    zplot = (
        np.power(yy, 2.)
        * ed.dtdfp(mc, xx, e_of_fp_interp(xx))
    )

    plt.figure(figsize=(15, 15))
    plt.contourf(xx, yy, zplot)
    plt.xscale('log')
    plt.colorbar()
    plt.xlabel(r'$f_p$')
    plt.ylabel(r'$\chi$ [Mpc]')
    plt.title(r'$\frac{dt[f_p, m_c]}{df_p} r^2, \quad m_c = $' + f'{mc:.2f}' + r' $M_\odot$')
    plt.savefig(os.path.join(savepath, 'r2dtdfp.png'))
    
    zplot = (
        ed.dtdfp(mc, xx, e_of_fp_interp(xx))
    )

    plt.figure(figsize=(15, 15))
    plt.contourf(xx, yy, zplot)
    plt.xscale('log')
    plt.colorbar()
    plt.xlabel(r'$f_p$')
    plt.ylabel(r'$\chi$ [Mpc]')
    plt.title(r'$\frac{dt[f_p, m_c]}{df_p}, \quad m_c = $' + f'{mc:.2f}' + r' $M_\odot$')
    plt.savefig(os.path.join(savepath, 'dtdfp.png'))










