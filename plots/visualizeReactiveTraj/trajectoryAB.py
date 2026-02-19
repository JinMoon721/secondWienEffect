import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as mpimg

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
lsize=19
mpl.rcParams['font.size']=fsize
#mpl.rcParams['pdf.fonttype']=42
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


rdata=load_ragged("results/LiPF6inACND05E45AB.r")
adata=load_ragged("results/LiPF6inACND05E45AB.a")

#rdata=load_ragged("results/ImpX20D010E50AB.r")
#adata=load_ragged("results/ImpX20D010E50AB.a")

colors=sns.color_palette("rocket", 3)

cc=sns.color_palette("deep", 20)
figA=1

fig, ax= plt.subplots(figsize=(10,5))

if figA:

    every=1

    cs = np.linspace(-1, 1, num=100);
    sn = np.sqrt(1-cs**2)
    r=5
    plt.plot(r*cs, r*sn, '--k')
    r=9
    plt.plot(r*cs, r*sn, '--k')
    
    count=0
    for tj in np.arange(1200, 5000, 267):
        x=rdata[tj] * adata[tj];
        y=rdata[tj] * np.sqrt( 1-adata[tj]**2)
        plt.plot(x, y, lw=3, color=cc[count])
        count+=1

    nx=[]
    ny=[]
    for tj in range(len(rdata)):
        nx.append(rdata[tj][0] * adata[tj][0])
        ny.append(rdata[tj][0] * np.sqrt(1-adata[tj][0]**2))
    plt.plot(nx, ny, 'o', color=colors[1], ms=0.5)

    nx=[]
    ny=[]
    for tj in range(len(rdata)):
        nx.append(rdata[tj][-1] * adata[tj][-1])
        ny.append(rdata[tj][-1] * np.sqrt(1-adata[tj][-1]**2))
    plt.plot(nx, ny, 'o', color=colors[0], ms=0.5)



ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect("equal", adjustable="box")
ax.axis('off')

img=mpimg.imread("../pictures/AB.png")
axImg = fig.add_axes([0.435, 0.1, 0.25, 0.25])
axImg.imshow(img)
axImg.axis("off")
    
plt.savefig('figure1.png', dpi=300, bbox_inches="tight")
plt.show()
