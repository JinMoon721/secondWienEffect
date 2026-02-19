# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 14:31:27 2025

@author: whisk
"""


from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE
import matplotlib.image as mpimg

#fig = plt.figure(1, figsize=(7.3, 4))
fig = plt.figure(1, figsize=(7.3, 4.))

plt.tick_params(direction='in', right=True, top=True)
plt.tick_params(labelsize=fsize)
plt.minorticks_on()
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)


colors=sns.color_palette("rocket", 5)
sigma=2

density="05"
target="LiPF6inACN"
filename="./results/rate/rate%sD%s" % (target, density)

pid=0
lb1=r"$\mathrm{CH_3CN}$"
mk="s"
df = pd.read_csv(filename, sep=r"\s+", header=None, comment="#", 
                 names=["density", "field", "cutin", "cutout", "q+", "Eq+", "q-", "Eq-", 
                        "flux", "Eflux", "kAB", "EkAB", "kBA", "EkBA", "kbm", "Ekbm",
                        "pA", "pB", "pAB", "mfptAB", "EmfptAB", "mfptBA", "EmfptBA"
                        ])

df["field"]*=convert ## from kcal/molA to mV/A
field=df["field"].to_numpy()[index]
dis=df["kAB"].to_numpy()[index]
disE=df["EkAB"].to_numpy()[index]

col="#D55E00"
y=np.log(dis/dis[0])
plt.plot(field, y, linestyle='-', lw=3, marker=mk, color=col, mfc='w', markersize=MS, label=lb1)
error = np.zeros((len(dis)))
for i in range(len(dis)):
    error[i] = divisionE(dis[i], disE[i], dis[0], disE[0])*sigma
plt.fill_between(field , y+error/np.exp(y), y-error/np.exp(y), alpha=AP, color=col)


target="LiPF6inH2O"
filename="./results/rate/rate%sD%s" % (target, density)

pid=1
lb1=r"$\mathrm{H_2O}$"
mk="o"
df = pd.read_csv(filename, sep=r"\s+", header=None, comment="#", 
                 names=["density", "field", "cutin", "cutout", "q+", "Eq+", "q-", "Eq-", 
                        "flux", "Eflux", "kAB", "EkAB", "kBA", "EkBA", "kbm", "Ekbm",
                        "pA", "pB", "pAB", "mfptAB", "EmfptAB", "mfptBA", "EmfptBA"
                        ])

df["field"]*=convert ## from kcal/molA to mV/A
field=df["field"].to_numpy()[index]
dis=df["kAB"].to_numpy()[index]
disE=df["EkAB"].to_numpy()[index]

col="#2A9D8F"
y=np.log(dis/dis[0])
plt.plot(field, y, linestyle='-', lw=3, marker=mk, color=col, mfc='w', markersize=MS, label=lb1)
error = np.zeros((len(dis)))
for i in range(len(dis)):
    error[i] = divisionE(dis[i], disE[i], dis[0], disE[0])*sigma
plt.fill_between(field , y+error/np.exp(y), y-error/np.exp(y), alpha=AP, color=col)

plt.plot(field, field*0, '--', lw=3, color='k')

plt.xlabel(r'$\xi \ / \ \mathrm{mV \ \AA^{-1}}$', fontsize=fsize)
plt.ylabel(r'$\ln k_{\mathrm{a}}(\xi) / k_{\mathrm{a}}(0) $', fontsize=fsize)
plt.xlim(0, 51)
plt.ylim(-0.2, 0.8)
plt.legend(fontsize=lsize, frameon=False)

img=mpimg.imread("pictures/association.png")
axImg = fig.add_axes([0.14, 0.48, 0.35, 0.35])
axImg.imshow(img)
axImg.axis("off")


plt.show()