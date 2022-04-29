import copy 
import eccentricdark as ed
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import sys

sns.set()

fig1 = False 
fig2 = False 
fig3 = False 
fig4 = False
extra1 = True
extra2 = False 

#Collect output directory 
if(len(sys.argv)!=2): 
    sys.exit("\nERROR: Please (only) specify the output directory when calling this script!\n") 
else: 
    savepath = sys.argv[1]
    print("\nSaving results to: ", savepath, "\n")

#Check LISA SNR - Figure 1 
if fig1==True:
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


# Check primordial distribution functions - Figure 2
if fig2==True: 
    nsample = 10000
    isolated_sample = ed.estar_sampler("isolated")(np.random.random(nsample))
    ejected_sample = ed.estar_sampler("ejected")(np.random.random(nsample))
    incluster_sample = ed.estar_sampler("incluster")(np.random.random(nsample))
    galcenter_sample = ed.estar_sampler("galcenter")(np.random.random(nsample))

    edges = np.logspace(-7, -2, 51)
    isolated_counts, bins = np.histogram(isolated_sample, edges)
    ejected_counts, bins = np.histogram(ejected_sample, edges)
    incluster_counts, bins = np.histogram(incluster_sample, edges)
    galcenter_counts, bins = np.histogram(galcenter_sample, edges)

    isolated_interp = scipy.interpolate.interp1d(edges[0:-1], isolated_counts/max(isolated_counts))
    ejected_interp = scipy.interpolate.interp1d(edges[0:-1], ejected_counts/max(ejected_counts))
    incluster_interp = scipy.interpolate.interp1d(edges[0:-1], incluster_counts/max(incluster_counts))
    galcenter_interp = scipy.interpolate.interp1d(edges[0:-1], galcenter_counts/max(galcenter_counts))

    plt.figure(figsize=(26.25, 7.5))
    plt.hist(bins[:-1], bins, weights=isolated_counts/max(isolated_counts), 
        label=r'isolated', linewidth=2, histtype='stepfilled', color='grey', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=ejected_counts/max(ejected_counts), 
        label=r'ejected', linewidth=2, histtype='stepfilled', color='red', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=incluster_counts/max(incluster_counts), 
        label=r'in-cluster', linewidth=2, histtype='stepfilled', color='purple', alpha=0.5)
    plt.hist(bins[:-1], bins, weights=galcenter_counts/max(galcenter_counts), 
        label=r'gal. center', linewidth=2, histtype='stepfilled', color='blue', alpha=0.5)
    plt.plot(np.power(10., ed.fieldData[:,0]), ed.fieldData[:,1]/max(ed.fieldData[:,1]), 
        linewidth=4, color="grey")
    plt.plot(np.power(10., ed.ejectedData[:,0]), ed.ejectedData[:,1]/max(ed.ejectedData[:,1]), 
        linewidth=4, color="red")
    plt.plot(np.power(10., ed.inclusterData[:,0]), ed.inclusterData[:,1]/max(ed.inclusterData[:,1]), 
        linewidth=4, color="purple")
    plt.plot(np.power(10., ed.galcenterData[:,0]), ed.galcenterData[:,1]/max(ed.galcenterData[:,1]), 
        linewidth=4, color="blue")
    plt.plot(np.logspace(-7, -4, 31), ed.multigauss(
            y = np.logspace(-7, -4, 31),
            u = -1.24404586e+01,
            a = -4.34697804e-01,
            b = -1.01499373e-01,
            c = -1.35644612e-02,
            k = 1.66563843e-06
        ), color="grey", linestyle='dashed', linewidth=4)
    plt.plot(np.logspace(-7, -3, 31), ed.multigauss(
            y = np.logspace(-7, -3, 31),
            u = -8.78221641e+00,
            a = -2.09640376e-01,
            b = -5.04274365e-02,
            c = -5.23793259e-03,
            k = 9.08204226e-06 
        ), color="red", linestyle='dashed', linewidth=4)
    plt.plot(np.logspace(-7, -2, 31), ed.multigauss(
            y = np.logspace(-7, -2, 31),
            u = -7.19800218e+00,
            a = -3.87340001e-01,
            b = -1.24883309e-01,
            c = -1.69861077e-02,
            k = 1.61038230e-04
        ), color="purple", linestyle='dashed', linewidth=4)
    plt.plot(np.logspace(-7, -2, 31), ed.multigauss(
            y = np.logspace(-7, -2, 31),
            u = -4.74123592e+00,
            a = -5.54338733e-01,
            b = -1.76646912e-01,
            c = -1.76665852e-02,
            k = 1.01316664e-03
        ), color="blue", linestyle='dashed', linewidth=4)

    plt.legend(fontsize=20)
    plt.ylabel(r"$e_0$-distribution", fontsize=20)
    plt.xlabel(r"$e_0$ at $f_{p0}=10$ Hz", fontsize=20)
    plt.xscale('log') 
    plt.xlim((10**-7, 10**-2))
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










