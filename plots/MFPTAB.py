# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 14:31:27 2025

@author: whisk
"""


from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE
import matplotlib.image as mpimg

fig = plt.figure(1, figsize=(5.5, 4.5))
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
disE=df["EkAB"].to_numpy()[index] * sigma

mfpt=df["mfptAB"].to_numpy()[index]
mfptE=df["EmfptAB"].to_numpy()[index] * sigma


mk="p"
lb1=r"$\nu_{AB}/(\langle q_- \rangle)$"
plt.plot(field, dis, linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=MS, label=lb1)
plt.fill_between(field , dis+disE, dis-disE, alpha=AP, color=colors[pid])

mk="D"
pid+=1
lb1=r"$\langle \tau_{AB} \rangle^{-1}$"
plt.plot(field, mfpt, linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=MS, label=lb1)
plt.fill_between(field , mfpt+mfptE, mfpt-mfptE, alpha=AP, color=colors[pid])


plt.xlabel(r'$\xi \ / \ \mathrm{mV \cdot \AA^{-1}}$', fontsize=fsize)
plt.ylabel(r'$k_{\mathrm{asso}}(\xi)  \ / \ \mathrm{ns^{-1}}$', fontsize=fsize)
plt.xlim(0, 51)
plt.legend(fontsize=lsize, frameon=False)

#img=mpimg.imread("pictures/dissociation.png")
#axImg = fig.add_axes([0.15, 0.38, 0.40, 0.40])
#axImg.imshow(img)
#axImg.axis("off")


plt.show()