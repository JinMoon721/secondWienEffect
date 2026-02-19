import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.tri as mtri

fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=15
mpl.rcParams['font.size']=fsize
mpl.rcParams['font.family']='serif'
mpl.rc('text', usetex=True)

import pandas as pd

if 1:
    fig, ax = plt.subplots(figsize=(7, 5.))
    gs=gridspec.GridSpec(100, 100)
    gs.update(wspace=5, hspace=2)
    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    colors=sns.color_palette("rocket", 10)

filenameK00 = "../processedData/opt/scanKPF6inACND05E00"
filenameK40 = "../processedData/opt/scanKPF6inACND05E40"
filenameLi00 = "../processedData/opt/scanLiPF6inACND05E00"
filenameLi40 = "../processedData/opt/scanLiPF6inACND05E40"

dfK00 = pd.read_csv(filenameK00, sep=r"\s+", header=None, comment="#")
dfK40 = pd.read_csv(filenameK40, sep=r"\s+", header=None, comment="#")
dfLi00 = pd.read_csv(filenameLi00, sep=r"\s+", header=None, comment="#")
dfLi40 = pd.read_csv(filenameLi40, sep=r"\s+", header=None, comment="#")

Li00cond = 66.527290
Li00condE = 2*0.934293

Li40cond = 90.012869
Li40condE = 0.85956995

K00cond = 71.699615
K00condE = 2*1.316145

K40cond = 94.736273
K40condE = 1.3112669

ratio00 = K00cond/Li00cond
ratio00E = ratio00 * ( (K00condE/K00cond)**2 + (Li00condE/Li00cond)**2)**(0.5)

ratioK = K40cond/K00cond
ratioKE = ratioK * ( (K40condE/K40cond)**2 + (K00condE/K00cond)**2)**(0.5)

ratioLi = Li40cond/Li00cond
ratioLiE = ratioLi * ( (Li40condE/Li40cond)**2 + (Li00condE/Li00cond)**2)**(0.5)

tol = 1e-3
rin=np.linspace(4.3, 5.9, num=9)

if 1:
    Liratio = dfLi40[6] / dfLi00[6]
    chiLi = (Liratio - ratioLi)/ratioLiE

    Kratio = dfK40[6] / dfK00[6]
    chiK = (Kratio - ratioK)/ratioKE

    ratio = dfK00[6] / dfLi00[6]
    chi00 = (ratio - ratio00)/ratio00E

    chisq = chiLi**2 + chiK**2 + chi00**2

    x=dfLi00[2].to_numpy() ## rin
    y=dfLi00[3].to_numpy() ## rout
    z=chisq.to_numpy()

    z=np.log10(z)

    mask = (z<2.5)

    x, y, z = x[mask], y[mask], z[mask]
    tri = mtri.Triangulation(x, y)
    interp = mtri.CubicTriInterpolator(tri, z)
    nx, ny = 300, 300
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    Xi, Yi = np.meshgrid(xi, yi)

    Zi = interp(Xi, Yi)

    cf = ax.contourf(Xi, Yi, Zi, levels=100)
    fig.colorbar(cf, ax=ax, label=r"$\log_{10} (\chi^2)$")

    plt.xlim(4.5, 5.8)
    plt.ylim(8, 11)
    plt.xlabel(r'$R_{\mathrm{i}} \ / \ \mathrm{\AA}$', fontsize=fsize)
    plt.ylabel(r'$R_{\mathrm{o}} \ / \ \mathrm{\AA}$', fontsize=fsize)
    

plt.savefig('./figures/sfigure1.png', dpi=300, bbox_inches="tight")
plt.show()
