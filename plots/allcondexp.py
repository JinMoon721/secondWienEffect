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
ndata=make_array(get_data("./results/conductivity/cond%sD%s" % (target, density)))
ahfield=ndata[:,1]
ahcond=ndata[:,4]
ahEcond=ndata[:,5]

target="LiPF6inH2O"
density="05"
ndata=make_array(get_data("./results/conductivity/cond%sD%s" % (target, density)))
hfield=ndata[:,1]
hcond=ndata[:,4]
hEcond=ndata[:,5]

target="KPF6inACN"
density="05"
ndata=make_array(get_data("./results/conductivity/cond%sD%s" % (target, density)))
kfield=ndata[:,1]
kcond=ndata[:,4]
kEcond=ndata[:,5]

target="LiPF6inACN"
density="05"
ndata=make_array(get_data("./results/rate/rate%sD%s" % (target, density)))
ahqn=ndata[:,6]
ahqnE=ndata[:,7]

target="LiPF6inH2O"
density="05"
ndata=make_array(get_data("./results/rate/rate%sD%s" % (target, density)))
hqn=ndata[:,6]
hqnE=ndata[:,7]

convert = 25.7/0.592
qx=ndata[:,1] * convert

target="KPF6inACN"
density="05"
ndata=make_array(get_data("./results/rate/rate%sD%s" % (target, density)))
kqn=ndata[:,6]
kqnE=ndata[:,7]

kx=ndata[:,1] * convert


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

    

colors=sns.color_palette("rocket", 3)

figA=1


if figA:
    fig, ax = plt.subplots(figsize=(8, 6))
    sigma=2

    index = list(range(22))
    skip = [2, 4, 6, 8, 10]
    index = [v for i, v in enumerate(index) if i not in skip]

    x=ahfield[index]
    y=ahcond[index]
    ey=ahEcond[index]
    qn=ahqn[index]
    qnE=ahqnE[index]
    qx=qx[index]
    ind=0
    ma='s'
    ms=8
    lb=r"$\lambda^\ast_\mathrm{ACN}$"
    lb2=r"$\langle q_- \rangle^\ast_\mathrm{ACN}$"
    field, cond, Econd, zerofield= diff(x, y, ey)

    plt.plot(field, 100*(cond-1), linestyle='-', lw=0,  marker=ma, color=colors[ind], mfc=colors[ind], markersize=ms, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=colors[ind], capsize=3, elinewidth=1, linewidth=1, ecolor=colors[ind], markeredgecolor=colors[ind])

    error = np.zeros((len(qn)-1))
    for i in range(len(qn)-1):
        error[i] = qn[i+1]/qn[0] * np.sqrt( (qnE[i+1]/qn[i+1])**2 + (qnE[0]/qn[0])**2)
    plt.plot(qx[1:], 100*(qn[1:]/qn[0]-1) , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label=lb2)
    plt.fill_between(qx[1:], 100*(qn[1:]/qn[0]-1 + error*sigma), 100*(qn[1:]/qn[0] - error*sigma-1), alpha=0.5, color=colors[ind])


    x=hfield[index]
    y=hcond[index]
    ey=hEcond[index]
    qn=hqn[index]
    qnE=hqnE[index]

    ind=1
    ma='o'
    lb=r"$\lambda^\ast_\mathrm{H_2O}$"
    lb2=r"$\langle q_- \rangle^\ast_\mathrm{H_2O}$"
    field, cond, Econd, zerofield= diff(x, y, ey)

    plt.plot(field, 100*(cond-1), linestyle='-', lw=0, marker=ma, color=colors[ind], mfc=colors[ind], markersize=ms, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=colors[ind], capsize=3, elinewidth=1, linewidth=1, ecolor=colors[ind], markeredgecolor=colors[ind])

    plt.plot(qx[1:], 100*(qn[1:]/qn[0]-1) , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label=lb2)
    error = np.zeros((len(qn)-1))
    for i in range(len(qn)-1):
        error[i] = qn[i+1]/qn[0] * np.sqrt( (qnE[i+1]/qn[i+1])**2 + (qnE[0]/qn[0])**2)
    plt.fill_between(qx[1:], 100*(qn[1:]/qn[0]-1 + error*sigma), 100*(qn[1:]/qn[0]-1 - error*sigma), alpha=0.5, color=colors[ind])


    index = list(range(21))
    skip = [2, 4, 6, 8, 10]
    index = [v for i, v in enumerate(index) if i not in skip]

    x=kfield[index]
    y=kcond[index]
    ey=kEcond[index]
    qn=kqn[index]
    qnE=kqnE[index]
    qx = kx[index]

    ind=2
    ma='p'
    lb=r"$\lambda^\ast_\mathrm{K}$"
    lb2=r"$\langle q_- \rangle^\ast_\mathrm{K}$"
    field, cond, Econd, zerofield= diff(x, y, ey)

    plt.plot(field, 100*(cond-1), linestyle='-', lw=0, marker=ma, color=colors[ind], mfc=colors[ind], markersize=ms, label=lb)
    ax.errorbar(field, 100*(cond-1), yerr = 100*Econd, fmt=ma, mfc=colors[ind], capsize=3, elinewidth=1, linewidth=1, ecolor=colors[ind], markeredgecolor=colors[ind])

    plt.plot(qx[1:], 100*(qn[1:]/qn[0]-1) , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label=lb2)
    error = np.zeros((len(qn)-1))
    for i in range(len(qn)-1):
        error[i] = qn[i+1]/qn[0] * np.sqrt( (qnE[i+1]/qn[i+1])**2 + (qnE[0]/qn[0])**2)
    plt.fill_between(qx[1:], 100*(qn[1:]/qn[0]-1 + error*sigma), 100*(qn[1:]/qn[0]-1 - error*sigma), alpha=0.5, color=colors[ind])

    plt.plot(field, field*0, 'k--')
#    plt.xlim(0, 50)
#    plt.ylim(-2,3.5)

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    plt.xlabel(r'$\xi \ / \ \mathrm{mV \cdot \AA^{-1}}$', fontsize=fsize)
    #plt.ylabel(r'$\lambda(\xi) \ / \ \mathrm{S \cdot cm^2 \cdot mol^{-1}} $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, ncol=2)
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
    
plt.savefig('./figures/cond%s.png' % (target), dpi=300, bbox_inches="tight")
plt.show()
