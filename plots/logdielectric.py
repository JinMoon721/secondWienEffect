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
lsize=22
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


#target="Stock"

density="05"
target="LiPF6inACN"
ndata=make_array(get_data("./results/dielectric/diel%sD%s" % (target, density)))
hfield=ndata[:,1]
hsave=ndata[0,2] ## from fluctuation
hdielec=ndata[:,4]
hEdielec=ndata[:,5]

target="LiPF6inH2O"
ndata=make_array(get_data("./results/dielectric/diel%sD%s" % (target, density)))
lfield=ndata[:,1]
lsave=ndata[0,2] ## from fluctuation
print(lsave)
ldielec=ndata[:,4]
lEdielec=ndata[:,5]

convert = 25.7/0.592


hdielec *= hfield * convert
hEdielec *= hfield * convert
ldielec *= lfield * convert
lEdielec *= lfield * convert

colors=sns.color_palette("rocket", 5)

figA=1

fig = plt.figure(1, figsize=(8,6))
gs=gridspec.GridSpec(100, 100)
gs.update(wspace=5, hspace=2)

if figA:
    xtr_subplot=fig.add_subplot(gs[0:100, 0:100])

    st=1
    ms=8
    dydx = np.gradient(hdielec, hfield*convert, edge_order=2)
    dydx, Edydx = gradient_with_uncertainty(hfield*convert, hdielec, hEdielec)
    Edydx*=3
    plt.plot(hfield[st:] * convert, dydx[st:]-1, linestyle='-', lw=3, marker='s', color=colors[0], mfc='w', markersize=ms, label=r"$\mathrm{ACN}$")
    plt.fill_between(hfield[st:] * convert, dydx[st:]+Edydx[st:]-1, dydx[st:]-Edydx[st:]-1, alpha=0.5, color=colors[0])
    plt.plot(hfield*convert, hfield*0 + hsave, '--', color=colors[0])

    dydx = np.gradient(ldielec, lfield*convert, edge_order=2)
    dydx, Edydx = gradient_with_uncertainty(lfield*convert, ldielec, lEdielec)
    Edydx*=3
    plt.plot(lfield[st:] * convert, dydx[st:]-1, linestyle='-', lw=3, marker='o', color=colors[1], mfc='w', markersize=ms, label=r"$\mathrm{H_2O}$")
    plt.fill_between(lfield[st:] * convert, dydx[st:]+Edydx[st:]-1, dydx[st:]-Edydx[st:]-1, alpha=0.5, color=colors[1])
    plt.plot(lfield*convert, lfield*0 + lsave, '--', color=colors[1])


    x=hfield[15:]*convert
    NA= 6.022 * 10**(23)
    rho = 55.5 * NA / 10**(-3)
    #rho = 15.5 * NA / 10**(-3)

    mu=2.0 * 3.335* 10**(-30) 
    epsilon0 = 8.854187* 10**(-12)
    beta = NA/0.592 / 4184
    xx = x * 10**(7)

#    plt.plot(x, 1+1/epsilon0*rho / beta /xx**2, '--k', label=r"$1+\frac{k_B T \rho}{\epsilon_0 \xi^2}$")

    plt.plot(x, 10**5/x**2, '--')
    x=hfield[1:25]*convert
    plt.plot(x, 10**3/x, '--')


    plt.xscale("log")
    plt.yscale("log")

#    plt.xlim(0, 50)
#    plt.ylim(-2,3.5)

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    plt.xlabel(r'$\xi \ / \ \mathrm{mV \cdot \AA^{-1}}$', fontsize=fsize)
    plt.ylabel(r'$\epsilon_r $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False)
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
    
plt.savefig('./figures/dielec.png' , dpi=300, bbox_inches="tight")
plt.show()
