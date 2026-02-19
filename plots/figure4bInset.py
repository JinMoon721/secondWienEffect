import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.image as image
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as mpimg



fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=19
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

rdata=load_ragged("../processedData/mechanism/ImpX20D010E50BA.r")
adata=load_ragged("../processedData/mechanism/ImpX20D010E50BA.a")

colors=sns.color_palette("rocket", 3)
cc=sns.color_palette("deep", 20)

fig, ax= plt.subplots(figsize=(10,5))
if 1:
    cs = np.linspace(-1, 1, num=1000);
    sn = np.sqrt(1-cs**2)
    r=5
    plt.plot(r*cs, r*sn, '--k')
    r=9
    plt.plot(r*cs, r*sn, '--k')
    
    count=0
    for tj in np.arange(1200, 5000,210 ):
        x=rdata[tj] * adata[tj];
        y=rdata[tj] * np.sqrt( 1-adata[tj]**2)
        plt.plot(x, y, lw=3, color=cc[count] )
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

img=mpimg.imread("./pictures/BA.png")
axImg = fig.add_axes([0.339, 0.1, 0.25, 0.25])
axImg.imshow(img)
axImg.axis("off")
    
plt.savefig('./figures/figure4bInset.png', dpi=300, bbox_inches="tight")
plt.show()
