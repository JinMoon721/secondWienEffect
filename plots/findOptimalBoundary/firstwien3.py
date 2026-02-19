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
ndata=make_array(get_data("./results/rateC%sD%s" % (target, density)))
Lqn=ndata[:,6]
LqnE=ndata[:,7]

target="KPF6inACN"
density="05"
ndata=make_array(get_data("./results/rateC%sD%s" % (target, density)))
Kqn=ndata[:,6]
KqnE=ndata[:,7]

cutoff=ndata[:,3]


convert = 25.7/0.592
    

colors=sns.color_palette("rocket", 5)

figA=1


if figA:
    fig, ax = plt.subplots(figsize=(8, 6))
    sigma=2


    ma="s"
    ind=1
    ms=8

#    plt.plot(cutoff, qna , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label="ACN")
#    plt.plot(cutoff, qnh , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label="H2O")
    plt.plot(cutoff, Kqn/Lqn , linestyle='-', marker=ma, color=colors[ind], mfc='w', markersize=ms, label="K/L")

    lambdaK=74.93093872
    lambdaL=69.24489594

    plt.plot(cutoff, cutoff*0 + lambdaK/lambdaL, '--')


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
    plt.legend(fontsize=lsize, frameon=False)
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
    
#plt.savefig('./figures/cond%s.png' % (target), dpi=300, bbox_inches="tight")
plt.show()
