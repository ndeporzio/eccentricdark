import copy
import eccentricdark as ed
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import sys

sns.set()

#Collect output directory 
if(len(sys.argv)!=2):
    sys.exit("\nERROR: Please (only) specify the output directory when calling this script!\n")
else:
    savepath = sys.argv[1]
    print("\nSaving results to: ", savepath, "\n")

scale=1.0

data1 = np.loadtxt(savepath+"/nocosmo_chi3000uniform_estar3p5.txt")
data2 = np.loadtxt(savepath+"/nocosmo_chi3000uniform_estar4p0.txt")
data3 = np.loadtxt(savepath+"/nocosmo_chi3000uniform_estar4p5.txt")
data4 = np.loadtxt(savepath+"/nocosmo_chi3000uniform_estar5p0.txt")
data5 = np.loadtxt(savepath+"/nocosmo_chi3000uniform_estar5p5.txt")

data = [data1, data2, data3, data4, data5]

#Visualize observable number counts
colors = [
    (193/255., 178/255., 240/255., 0.7), #purple  
    (178/255., 209/255., 240/255., 0.7), #blue 
    (178/255., 240/255., 209/255., 0.7), #green 
    (255/255., 224/255., 178/255., 0.7) #orange
]

estar = np.array([-3.5, -4.0, -4.5, -5.0, -5.5])
e_cutoffs=np.array([0.01, 0.1, 0.4, 0.9])
fpbins = np.logspace(-2.7, -1.5, 13)
ecutorder = np.array([3,2,1,0])


for estaridx, estarval in enumerate(estar): 
    plt.figure(figsize=(7.5, 7.5))

    for ecutidx in ecutorder: 
        plt.hist(
            np.log10(fpbins)[0:-1],
            weights=scale*data[estaridx][ecutidx], 
            bins=np.log10(fpbins),
            color=colors[ecutidx],  
            label=r'$e < $'+f'{e_cutoffs[ecutidx]:.3f}',
            linewidth=4, 
            histtype='stepfilled'
        )

    plt.title("Fixed $e_*=$"+f'{estarval:.2e}'+", Static Cosmology", fontsize=20)
    plt.legend(fontsize=20)
    plt.xlabel(r"$\log (f_p/Hz)$", fontsize=20)
    plt.ylabel(r"LISA N2A5 Observable Counts", fontsize=20)
    plt.savefig(savepath+"/nocosmo_chivaried_estarfixed_"+f'{estarval:.2e}'+".png")
