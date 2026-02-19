import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

##
#import system as sys

#if (len(sys.argv) < 4):
#    print("Usage: python rdf.py target density")
#    sys.exit(1)

#target = sys.argv[1]
#print("Target : ", target)
#density = sys.argv[2]
#print("Density : ", density)




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
field="00"
ndata=make_array(get_data("./results/rdf/%sD%sE%s.dat" % (target, density, field)))
ar=ndata[:,0]
ardf=ndata[:,1]

target="LiPF6inH2O"
ndata=make_array(get_data("./results/rdf/%sD%sE%s.dat" % (target, density, field)))
hr=ndata[:,0]
hrdf=ndata[:,1]



colors=sns.color_palette("rocket", 3)


target="LiPF6inH2O"
density="05"
field="19"
ndata=make_array(get_data("./results/rdf/%sD%sE%s.dat" % (target, density, field)))
x1=ndata[:,0]
y1=ndata[:,1]


figA=1

fig = plt.figure(1, figsize=(6,5))
gs=gridspec.GridSpec(100, 100)
gs.update(wspace=5, hspace=2)

if figA:
    xtr_subplot=fig.add_subplot(gs[0:100, 0:100])

    ay = -np.log(ardf)
    aloc = ay.argmin()
    ay -= ay[aloc]

    y1 = -np.log(y1)
    yloc = y1.argmin()
    y1 -= y1[yloc]

    hy = -np.log(hrdf)
    hloc = hy.argmin()
    hy -= hy[hloc]

    kbT=0.592

    mask = (ar > 5.3)
    ay[mask] = ay[mask]-0.05/kbT
    plt.plot(ar,ay*kbT, linestyle='-',lw=3, marker='s', color=colors[0], mfc='w', markersize=0, label=r"$\mathrm{CH_3CN}$")
    plt.plot(x1,y1*kbT, linestyle='-',lw=3, marker='s', color=colors[2], mfc='w', markersize=0, label=r"$\mathrm{CH_3CN}$")


    mask = (hr > 5.3)
    hy[mask] = hy[mask]-0.05/kbT
    plt.plot(hr,hy*kbT, linestyle='-',lw=3, marker='o', color=colors[1], mfc='w', markersize=0, label=r"$\mathrm{H_2O}$")

    py = np.linspace(-0.5, 6, num=100)

    plt.plot(py*0 + 5.3, py, 'k--')
    plt.plot(py*0 + 9.1, py, 'k--')


    plt.xlim(2.4, 11.5)
    plt.ylim(-0.5,6)

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    plt.xlabel(r'$r \ / \ \mathrm{\AA}$', fontsize=fsize)
    plt.ylabel(r'$-k_B T \ln \rho(r) \ / \ \mathrm{kcal \cdot mol^{-1}} $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, loc='center')
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
#    plt.title(r"$\mathrm{LiPF_6 \ 0.5 \ M}$")
    
plt.savefig('./figures/rdf.png' , dpi=300, bbox_inches="tight")
plt.show()
