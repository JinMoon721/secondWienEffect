import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

def gradient_with_uncertainty(x, y, sigma_y):
    x = np.asarray(x); y = np.asarray(y); s = np.asarray(sigma_y)
    n = len(x)
    g  = np.empty(n)
    sg = np.empty(n)

    # interior points: 3-point central, non-uniform
    for i in range(1, n-1):
        h1 = x[i]   - x[i-1]
        h2 = x[i+1] - x[i]
        a_m1 = -h2/(h1*(h1+h2))
        a_0  =  (h2-h1)/(h1*h2)
        a_p1 =  h1/(h2*(h1+h2))
        g[i]  = a_m1*y[i-1] + a_0*y[i] + a_p1*y[i+1]
        sg[i] = np.sqrt((a_m1*s[i-1])**2 + (a_0*s[i])**2 + (a_p1*s[i+1])**2)

    # left edge: forward difference
    h = x[1] - x[0]
    g[0]  = (y[1]-y[0])/h
    sg[0] = np.sqrt(s[1]**2 + s[0]**2)/abs(h)

    # right edge: backward difference
    h = x[-1] - x[-2]
    g[-1]  = (y[-1]-y[-2])/h
    sg[-1] = np.sqrt(s[-1]**2 + s[-2]**2)/abs(h)

    return g, sg



fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=18
mpl.rcParams['font.size']=fsize
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['font.family']='serif'
mpl.rc('text', usetex=True)

def get_data(file_name):
    data=[]
    with open(file_name, "r") as f:
        lines=f.readlines()
        for i in range(len(lines)):
            l=[]
            for t in lines[i].split():
                try:
                    l.append(float(t))
                except ValueError:
                    pass
            data.append(l)
        data=np.asarray(data)
        data=data.astype(float)
    return data
def make_array(array):
    array=np.asarray(array)
    array=array.astype(float)
    return array


target="LiPF6inACN"
density="05"
ndata=make_array(get_data("../processedData/conductivity/cond%sD%s" % (target, density)))
ahfield=ndata[:,1]
ahcond=ndata[:,4]
ahEcond=ndata[:,5]


target="LiPF6inH2O"
density="05"
ndata=make_array(get_data("../processedData/conductivity/cond%sD%s" % (target, density)))
hfield=ndata[:,1]
hcond=ndata[:,4]
hEcond=ndata[:,5]

out="72" ## single cutoff at 7.2A
target="LiPF6inACN"
density="05"
ndata=make_array(get_data("../processedData/population/rout%s" % (out)))
pA=ndata[:,16]
pB=ndata[:,17]
pAB=ndata[:,18]

pop = pA
print(pA+pB+pAB)

convert = 25.7/0.592
qx=ndata[:,1] * convert
def diff(x, y, ey):
    zerofield=y[0]
    e2 = ey[0]/y[0]
    convert = 25.7/0.592
    sigma=2
    y*= convert*x
    ey *= convert*x
    dydx = np.gradient(y, x*convert, edge_order=1)
    dydx, Edydx = gradient_with_uncertainty(x*convert, y, ey)

    e1 = (Edydx[1:] / dydx[1:])
    error = dydx[1:]/zerofield * np.sqrt( e1**2 + e2**2)
    return x[1:]*convert, dydx[1:]/zerofield, sigma*error, zerofield

    
colors=sns.color_palette("rocket", 5)

if 1:
    fig, ax = plt.subplots(figsize=(7, 5))
    sigma=2

    ap=0.6
    x=ahfield
    y=ahcond
    ey=ahEcond
    pop=pop
    qx=qx

    ind=0
    ma='o'
    ms=10
    lb=r"$\lambda^\ast$"
    lb2=r"$p_A^\ast$"
    field, cond, Econd, zerofield= diff(x, y, ey)

    plt.plot(field, 100*(cond-1), linestyle='-', lw=0,  marker=ma, color=colors[ind], mfc=colors[ind], markersize=ms, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=colors[ind], capsize=3, elinewidth=1, linewidth=1, ecolor=colors[ind], markeredgecolor=colors[ind])

    plt.plot(qx[1:], 100*(pop[1:]/pop[0]-1) , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label=lb2)

    plt.plot(field, field*0, 'k--')
    plt.xlim(0, 53)
    plt.ylim(-17,52)

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)


    plt.xlabel(r'$\xi \ / \ \mathrm{mV \ \AA^{-1}}$', fontsize=fsize)
    plt.ylabel(r'$x^* = 100 \times (x_\xi - x_0)/x_0  $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, ncol=2, columnspacing=0.6)
    
plt.savefig('./figures/sfigure2a.png',  dpi=300, bbox_inches="tight")
plt.show()
