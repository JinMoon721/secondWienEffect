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


colacn="#D55E00"
colh2o="#2A9D8F"
fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=19
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

target="LiPF6inACN"
density="05"
ndata=make_array(get_data("../processedData/rate/rate%sD%s" % (target, density)))
ahqn=ndata[:,6]
ahqnE=ndata[:,7]

target="LiPF6inH2O"
density="05"
ndata=make_array(get_data("../processedData/rate/rate%sD%s" % (target, density)))
hqn=ndata[:,6]
hqnE=ndata[:,7]



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
#    return x[1:]*convert, dydx[1:], sigma*Edydx[1:], zerofield

    

colors=sns.color_palette("rocket", 5)

figA=1


if figA:
    fig, ax = plt.subplots(figsize=(7.3, 5.))
    sigma=2

    index = list(range(22))
    skip = [2, 4, 6, 8]
    index = [v for i, v in enumerate(index) if i not in skip]
    
    ap=0.6
    x=ahfield[index]
    y=ahcond[index]
    ey=ahEcond[index]
    qn=ahqn[index]
    qnE=ahqnE[index]
    qx=qx[index]
    ind=0
    ma='s'
    ms=10
    lb=r"$\lambda^\ast_\mathrm{CH_3CN}$"
    lb2=r"$\langle q_- \rangle^\ast_\mathrm{CH_3CN}$"
    field, cond, Econd, zerofield= diff(x, y, ey)
    col=colacn

    plt.plot(field, 100*(cond-1), linestyle='-', lw=0,  marker=ma, color=col, mfc=col, markersize=ms)#, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=col, capsize=5, elinewidth=2, linewidth=1, ecolor=col, markeredgecolor=col)

    error = np.zeros((len(qn)-1))
    for i in range(len(qn)-1):
        error[i] = qn[i+1]/qn[0] * np.sqrt( (qnE[i+1]/qn[i+1])**2 + (qnE[0]/qn[0])**2)
    plt.plot(qx[1:], 100*(qn[1:]/qn[0]-1) , linestyle='-', marker='o', color=col, mfc='w', markersize=ms)#, label=lb2)
    plt.fill_between(qx[1:], 100*(qn[1:]/qn[0]-1 + error*sigma), 100*(qn[1:]/qn[0] - error*sigma-1), alpha=ap, color=col)



    ### onsager's expression
    field = np.linspace(0.0001, 1.7, num=500)
    dielectric=20
    kBT=0.592
    qqrd2e =332.06371 ## e^2/4 pi epsilon0 in kcal/mol A unit
    dielectric=19.202
    kcal2mev = 25.7/kBT
    qE = field/kcal2mev ## kcal/molA 
    elength=kBT/qE ## A unit
    Bjerrum = qqrd2e/kBT/dielectric ## A unit
    x=Bjerrum/elength/2 ## factor of 2 from Onsager's definition

    from scipy.special import i1
    ratio = i1(np.sqrt(8*x)) / np.sqrt(2*x)
    conc=0.5
    bq=0.401
    k0=bq**2 * conc / (1-bq)
    k=k0 *ratio

    newratio = (np.sqrt( k**2 + 4*k*conc) - k) / (np.sqrt(k0**2 + 4*k0*conc) - k0)

    plt.plot(field, (newratio-1)*100,'--', lw=3, color="#D55E00")#, label="Onsager")

    x=hfield[index]
    y=hcond[index]
    ey=hEcond[index]
    qn=hqn[index]
    qnE=hqnE[index]

    ind=1
    ma='s'
    lb=r"$\lambda^\ast_\mathrm{H_2O}$"
    lb2=r"$\langle q_- \rangle^\ast_\mathrm{H_2O}$"
    field, cond, Econd, zerofield= diff(x, y, ey)
    col=colh2o
    plt.plot(field, 100*(cond-1), linestyle='-', lw=0, marker=ma, color=col, mfc=col, markersize=ms)#, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=col, capsize=5, elinewidth=2, linewidth=1, ecolor=col, markeredgecolor=col)

    plt.plot(qx[1:], 100*(qn[1:]/qn[0]-1) , linestyle='-', marker='o', color=col, mfc='w', markersize=ms)#, label=lb2)
    error = np.zeros((len(qn)-1))
    for i in range(len(qn)-1):
        error[i] = qn[i+1]/qn[0] * np.sqrt( (qnE[i+1]/qn[i+1])**2 + (qnE[0]/qn[0])**2)
    plt.fill_between(qx[1:], 100*(qn[1:]/qn[0]-1 + error*sigma), 100*(qn[1:]/qn[0]-1 - error*sigma), alpha=ap, color=col)


    #plt.plot(field, field*0, 'k', lw=1.5)
    
    field = np.linspace(0.0001, 9, num=500)
    qE = field/kcal2mev
    elength=kBT/qE
    dielectric=69.564
    Bjerrum = qqrd2e/kBT/dielectric ## A unit
    x=Bjerrum/elength/2 ## factor of 2 from Onsager's definition

    ratio = i1(np.sqrt(8*x)) / np.sqrt(2*x)
    conc=0.5
    bq=0.585
    k0=bq**2 * conc / (1-bq)
    k=k0 *ratio

    newratio = (np.sqrt( k**2 + 4*k*conc) - k) / (np.sqrt(k0**2 + 4*k0*conc) - k0)
    plt.plot(field, (newratio-1)*100,'--', lw=3, color="#2A9D8F")# label='Onsager')
    
    
    plt.xlim(0, 51)
 #   plt.ylim(-2,3.5)

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    plt.xlabel(r'$\xi \ / \ \mathrm{mV \cdot \AA^{-1}}$', fontsize=fsize)
    plt.ylabel(r'$x^* = 100 \times (x_\xi - x_0)/x_0  $', fontsize=fsize)

#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)

x=np.linspace(100, 101, 2)
y=x
plt.plot(x, y, lw=0, marker='s', mfc='0.5', markersize=ms, color='0.5', label=r"$\lambda^*$")
plt.plot(x, y, lw=3, marker='o', mfc='w', markersize=ms, color='0.5', label=r"$\langle q_- \rangle^*$")
plt.plot(x, y, '--', lw=3, color='0.5', label=r"$\mathrm{Onsager}$")


plt.ylim(-15.5,50)
plt.legend(fontsize=lsize, frameon=False, loc='upper left', ncol=2, columnspacing=0.6)


#plt.savefig('./figures/cond.png',  dpi=300, bbox_inches="tight")
plt.show()
