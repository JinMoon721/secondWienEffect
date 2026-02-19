import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)


class histogram():
    def __init__(self, limits, binwidth):
        self.limits = limits
        self.binwidth = binwidth
        self.num_bins = int( (limits[1]-limits[0] ) / binwidth )
        self.vals = np.arange( limits[0], limits[1], binwidth) + binwidth /2
        self.histo = np.zeros(self.num_bins)
        self.num_samples =0
    def add_data(self, data):
        self.num_samples +=1
        if data >= self.limits[0] and data < self.limits[1]:
            bin = int( (data-self.limits[0] ) / self.binwidth )
            self.histo[bin] +=1
    def normalize(self):
        self.prob = self.histo / (self.num_samples * self.binwidth ) 

    def lineplot(self, name):
        plt.plot(self.vals, self.prob, 'o-', label=name)
        plt.ylabel('pdf', fontsize=14)
        plt.xlabel('variable', fontsize=14)

    def barplot(self):
        plt.bar(self.vals, self.prob, width=0.9 * self.binwidth, edgecolor='k', color='Orange')
        plt.ylabel('tdf', fontsize=14)
        plt.xlabel('variable', fontsize=14)


fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=22
mpl.rcParams['font.size']=fsize
mpl.rcParams['font.family']='serif'
mpl.rc('text', usetex=True)

def load_ragged(filename):
    rows =[]
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            arr = np.fromstring(line, sep=" ")
            rows.append(arr)
    return rows 


rdata=load_ragged("../processedData/mechanism/LiPF6inACND05E50BA.r")
adata=load_ragged("../processedData/mechanism/LiPF6inACND05E50BA.a")


colors=sns.color_palette("rocket", 3)

if 1:
    fig, ax = plt.subplots( figsize=(7.3,4.))
    xmax = 1.0
    xmin = -1.0
    histBin=(xmax - xmin)/20
    histA = histogram([xmin, xmax], histBin)
    histB = histogram([xmin, xmax], histBin)

    for tj in range(len(rdata)):
        histA.add_data( adata[tj][-1])
        histB.add_data( adata[tj][0])

    histA.normalize()
    mk=">"
    pid=0
    MS=10
    lb1=r"$R_{out}$"
    plt.plot(histA.vals, histA.prob, linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=MS, label=lb1)

    histB.normalize()
    mk="<"
    pid=1
    MS=10
    lb1=r"$R_{in}$"
    plt.plot(histB.vals, histB.prob, linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=MS, label=lb1)

    plt.ylim(0, 3 )

    plt.minorticks_on()
    plt.tick_params(direction='in', right=True, top=True)
    plt.tick_params(labelsize=nsize)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
    plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)

    plt.xlabel(r'$\cos \theta$', fontsize=fsize)
    plt.ylabel(r'$p (\cos \theta)$', fontsize=fsize)
    
plt.savefig('./figures/figure4a.png', dpi=300, bbox_inches="tight")
plt.show()
