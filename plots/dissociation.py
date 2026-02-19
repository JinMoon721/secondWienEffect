# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 14:31:27 2025

@author: whisk
"""


from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE
import matplotlib.image as mpimg

#fig = plt.figure(1, figsize=(7.3, 4)) ## original 7.5 5.5
fig = plt.figure(1, figsize=(7.3, 4)) ## original 7.5 5.5

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
dis=df["kBA"].to_numpy()[index]
disE=df["EkBA"].to_numpy()[index]

#dis=df["flux"].to_numpy()[index]
5#disE=df["Eflux"].to_numpy()[index]


col="#D55E00"

#plt.plot(field, dis/dis[0], linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=MS, label=lb1)
plt.plot(field, np.log(dis/dis[0]), linestyle='-', lw=3, marker=mk, color=col, mfc='w', markersize=MS, label=lb1)

error = np.zeros((len(dis)))
for i in range(len(dis)):
    error[i] = divisionE(dis[i], disE[i], dis[0], disE[0])*sigma
#plt.fill_between(field , dis/dis[0]+error, dis/dis[0]-error, alpha=AP, color=colors[pid])
y=np.log(dis/dis[0])
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
dis=df["kBA"].to_numpy()[index]
disE=df["EkBA"].to_numpy()[index]

#dis=df["flux"].to_numpy()[index]
#disE=df["Eflux"].to_numpy()[index]
col="#2A9D8F"
y=np.log(dis/dis[0])
plt.plot(field, y, linestyle='-', lw=3, marker=mk, color=col, mfc='w', markersize=MS, label=lb1)
error = np.zeros((len(dis)))
for i in range(len(dis)):
    error[i] = divisionE(dis[i], disE[i], dis[0], disE[0])*sigma
plt.fill_between(field , y+error/np.exp(y), y-error/np.exp(y), alpha=AP, color=col)

plt.plot(field, field*0, '--', lw=3, color='k')


plt.xlabel(r'$\xi \ / \ \mathrm{mV \cdot \AA^{-1}}$', fontsize=fsize)
plt.ylabel(r'$\ln k_{\mathrm{d}}(\xi) / k_{\mathrm{d}}(0) $', fontsize=fsize)
plt.xlim(0, 51)
plt.ylim(-0.2, 0.8)
plt.legend(fontsize=lsize, loc='center right', frameon=False)

#plt.yscale("log")
a=0.03
#plt.plot(field, np.log(np.cosh(a*field)), '--', lw=5)

plt.plot(field, field/60-0.12, '--', color="#D55E00", lw=2)
plt.plot(field, field/190-0.085, '--', color=col, lw=2)


### onsager's expression
field = np.linspace(0.0001, 50, num=500)
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

#plt.plot(field, np.log(newratio),'--', lw=3, color="#D55E00")
#print(np.log(newratio[-1]))
dielectric=69.564
Bjerrum = qqrd2e/kBT/dielectric ## A unit
x=Bjerrum/elength/2 ## factor of 2 from Onsager's definition

ratio = i1(np.sqrt(8*x)) / np.sqrt(2*x)
conc=0.5
bq=0.585
k0=bq**2 * conc / (1-bq)
k=k0 *ratio

newratio = (np.sqrt( k**2 + 4*k*conc) - k) / (np.sqrt(k0**2 + 4*k0*conc) - k0)
#plt.plot(field, np.log(newratio),'--', lw=3, color="#2A9D8F")





#plt.plot(field, np.exp(field/190)/1.1, '--', lw=3)
#plt.plot(field, 1+field/60, lw=3)
#plt.yscale("log")
img=mpimg.imread("pictures/dissociation.png")
#axImg = fig.add_axes([0.26, 0.49, 0.25, 0.25])
axImg = fig.add_axes([0.14, 0.48, 0.35, 0.35])
axImg.imshow(img)
axImg.axis("off")


plt.show()