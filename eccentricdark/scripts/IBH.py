import copy
import eccentricdark as ed
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import sys
from matplotlib.patches import Rectangle

sns.set()

fig1 = False
fig2 = True

#Collect output directory 
if(len(sys.argv)!=2):
    sys.exit("\nERROR: Please (only) specify the output directory when calling this script!\n")
else:
    savepath = sys.argv[1]
    print("\nSaving results to: ", savepath, "\n")

# Check LISA SNR Functions
if fig1==True:
   # x = np.logspace(-5, 0, 501)
   # y = np.array([np.sqrt(ed.SnLISAdefault(xval)) for xval in x])
   # plt.figure(figsize=(15, 7.5))
   # plt.plot(x, y)
   # plt.grid(True, which='both', axis='both')
   # plt.xscale('log')
   # plt.yscale('log')
   # plt.xlabel(r"$f_p$ [Hz]", fontsize=20)
   # plt.ylabel(r"$(Noise)^{1/2}$", fontsize=20)
   # plt.xticks(fontsize=18)
   # plt.yticks(fontsize=18)
   # plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
   # plt.savefig(os.path.join(savepath, 'LISA_N2A5_1.png'))


   # x = np.logspace(-5, 0, 501)
   # y = np.array([np.power(xval, 0.5)*np.sqrt(ed.SnLISAdefault(xval)) for xval in x])
   # plt.figure(figsize=(15, 7.5))
   # plt.plot(x, y)
   # plt.grid(True, which='both', axis='both')
   # plt.xscale('log')
   # plt.yscale('log')
   # plt.xlabel(r"$f_p$ [Hz]", fontsize=20)
   # plt.ylabel(r"$(f_p \times Noise)^{1/2}$", fontsize=20)
   # plt.xticks(fontsize=18)
   # plt.yticks(fontsize=18)
   # plt.title("LISA N2A5 Configuration Noise Curve", fontsize=20)
   # plt.savefig(os.path.join(savepath, 'LISA_N2A5_2.png'))

    e0 = 0.
    obstime = 10.

    fp = np.logspace(-4, 1, 51)
    m = np.geomspace(1.0e1*ed.msun_in_kg, 1.0e6*ed.msun_in_kg, 6)
    y = np.array([[
            ed.roffmSNR8(100., f, e0, mval/2., mval/2., obstime)
            for f in fp]
            for mval in m]
        )
    y2 = np.nan_to_num(y)
    roffmSNR8intpl = scipy.interpolate.interp2d(fp, m, y2, kind='linear')

    evolutions = [[ed.BBHevolve(f, e0, mval/2., mval/2., obstime)
        for f in fp]
        for mval in m]

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
    #plt.xscale('log')
    #plt.yscale('log')
    plt.text(-2.5, -2.,"LISA")
    plt.text(-1.5, -2.,"TianQin")
    plt.text(-0.5, -2.,"DECIGO")
    plt.text(0.6, -2.,"LIGO")
    plt.ylim((-3., 6.))
    plt.xlim((-4, 1.))
    plt.grid(True, which='both', axis='both')
    plt.xlabel(r"log($f_p/$[Hz])", fontsize=20)
    plt.ylabel(r"log(r/[Mpc]) s.t. SNR>8", fontsize=20)
    plt.title("LISA N2A5 Sensitivity, $m_1 = m_2$", fontsize=20)
    plt.legend()
    plt.savefig(os.path.join(savepath, 'LISA_N2A5_3.png'))


#Check LISA SNR
if fig2==True:
    m1 = 5.0e4*ed.msun_in_kg
    m2 = 5.0e4*ed.msun_in_kg

    sns.set()
    fp_table = np.logspace(-4, 2, 61)
    estar_table = np.array([1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6])

    e_of_fp = np.array([[
        ed.e_solver(
            f_p=fpval,
            fp_star=10., #Units: Hz 
            e_star=estarval
        )
        for fpidx, fpval in enumerate(fp_table)]
        for estaridx, estarval in enumerate(estar_table)])

    plt.figure(figsize=(15, 7.5))
    plt.title("Fig. 1 in 1907.02283, m1=m2="+f'{m1/ed.msun_in_kg:.1e}'+r'$M_\odot$', fontsize=20)
    plt.xlabel(r"log($f_p$/[Hz])", fontsize=20)
    plt.ylabel(r"log($e$)", fontsize=20)

    palette = sns.color_palette("flare", len(estar_table))
    for estaridx, estarval in enumerate(estar_table):
        plt.plot(
            np.log10(fp_table), 
            np.log10(e_of_fp[estaridx, :]), 
            label=r"log($e_*$)="+f'{np.log10(estarval):.1f}', 
            color = palette[estaridx] 
        )
    plt.plot([-4, 2], [0., 0.], color='black', label=r"e=1.0")
    plt.plot([-4, 2], [-2.699, -2.699], color='black', alpha=0.2)

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
        1.523, 
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
    plt.text(1.5, -2.,"LIGO")


    t_coals = np.array([
        ed.year_in_seconds/12., 
        ed.year_in_seconds, 
        10.*ed.year_in_seconds,
        100.*ed.year_in_seconds
    ])
    labels=['1m', '1yr', '10yr', '100yr']
    e_stars = np.logspace(-8, -1, 71)
    fp_coals = np.zeros((len(t_coals), len(e_stars)))
    e_coals = np.zeros((len(t_coals), len(e_stars)))
    for t_idx, t_val in enumerate(t_coals): 
        print("Solving coalescence ", t_idx)
        for e_idx, e_val in enumerate(e_stars):
            evolution = ed.e_to_fp_interpolator(10., e_stars[e_idx])
            fp_of_e = evolution[2]
            e_of_fp = evolution[3]
            try:
                #fp = scipy.optimize.bisect(
                #    lambda f: ed.lifetime30Msun(f, e_of_fp(f)) - t_val,  
                #    evolution[0], 
                #    evolution[1]
                #) 
                fp = scipy.optimize.bisect(
                    lambda f: ed.tmerge(f, e_of_fp(f), m1, m2) - t_val,  
                    evolution[0], 
                    evolution[1]
                ) 
                e = e_of_fp(fp) 
                
                fp_coals[t_idx, e_idx] = fp
                e_coals[t_idx, e_idx] = e 
            except: 
                fp_coals[t_idx, e_idx] = 1.
                e_coals[t_idx, e_idx] = 1. 
                           
 
    for t_idx, t_val in enumerate(t_coals):
        plt.plot(
            np.log10(fp_coals[t_idx,:]), 
            np.log10(e_coals[t_idx,:]), 
            label=labels[t_idx], 
            color='magenta', 
            alpha = (t_idx+1)/len(t_coals)
        )      
        

    plt.legend(fontsize=20)
    plt.xlim((-4, 2))
    plt.ylim((-3, 0.1))
    plt.grid(True, which='both', axis='both')
    plt.savefig(os.path.join(savepath, 'Figure_1.png'))
