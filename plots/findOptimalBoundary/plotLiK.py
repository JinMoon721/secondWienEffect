import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=15
mpl.rcParams['font.size']=fsize
#mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['font.family']='serif'
mpl.rc('text', usetex=True)

import pandas as pd

if 1:
    fig= plt.figure(1, figsize=(8,6))
    gs=gridspec.GridSpec(100, 100)
    gs.update(wspace=5, hspace=2)
    xtr_subplot=fig.add_subplot(gs[0:100, 0:100])
    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    colors=sns.color_palette("rocket", 10)



filenameK = "./results/scanKPF6inACND05E00"
filenameLi = "./results/scanLiPF6inACND05E00"
dfK = pd.read_csv(filenameK, sep=r"\s+", header=None, comment="#")
dfLi = pd.read_csv(filenameLi, sep=r"\s+", header=None, comment="#")

Licond = 66.527290
LicondE = 2*0.934293

Kcond = 71.699615
KcondE = 2*1.316145

ratio = Kcond/Licond
ratioE = ratio * ( (KcondE/Kcond)**2 + (LicondE/Licond)**2)**(0.5)

tol = 1e-3
rin=np.linspace(4.5, 5.9, num=8)

for i in range(len(rin)):
    dfK_fixed = dfK[np.isclose(dfK[2], rin[i], atol=tol, rtol=0.0)] ## absolute tolerance
    dfLi_fixed = dfLi[np.isclose(dfLi[2], rin[i], atol=tol, rtol=0.0)] ## absolute tolerance
    x = dfLi_fixed[3] ## rout
    y = dfK_fixed[6] / dfLi_fixed[6] ## q- ratio
    ey = y * ( (dfLi_fixed[7]/dfLi_fixed[6])**2 + (dfK_fixed[7]/dfK_fixed[6])**2 )**(0.5)   ## error q-
    y=100*(y-1)
    plt.plot(x, y, linestyle='-', marker='s', color=colors[i], mfc='w', markersize=2, label=r"$r_{in}=%.1f \ \mathrm{\AA}$" % (rin[i]))
#    plt.fill_between(x, y-E00ey, y+ey, alpha=0.5,  color=colors[i])
    plt.plot(x, x*0+ 100*(ratio-1), '--', color=colors[0])
    plt.fill_between(x, 100*(ratio-1-ratioE), 100*(ratio-1+ratioE), alpha=0.3, color=colors[-1])


    
plt.xlabel(r'$r_{out} \ / \ \mathrm{\AA}$', fontsize=fsize)
plt.ylabel(r'$100 \ \times \frac{\Delta \langle q_- \rangle_\xi}{ \langle q_- \rangle_0}$', fontsize=fsize)
plt.legend(fontsize=lsize, frameon=False)
plt.savefig('figures/boundaryLiK.png', dpi=300, bbox_inches="tight")
plt.show()
