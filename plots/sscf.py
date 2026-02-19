# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 14:31:27 2025

@author: whisk
"""


from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE
import matplotlib.image as mpimg

fig = plt.figure(1, figsize=(7., 5.))
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
filename="./results/single/%sD%sE00.dat" % (target, density)

pid=0
lb1=r"$\mathrm{CH_3CN}$"
mk="s"
df = pd.read_csv(filename, sep=r"\s+", header=None, comment="#", 
                 names=["time", "sscf"
                        ])

x=df["time"].to_numpy()/1000
y=df["sscf"].to_numpy() / 0.3759


plt.plot(x, y, linestyle='-', lw=3, marker=mk, color=colors[pid], mfc='w', markersize=0, label=r"$C_{AB}$")
xx=np.logspace(-4.4, -3.3, num=100)
plt.plot(xx, xx*10**(3.3), '--', lw=3, color=colors[1], label=r"$C_{AB} \propto t^{1.0}$")
yy=np.logspace(-4.4, -1., num=100)
plt.plot(yy, yy**0.4 *10**(.3), '--', lw=3, color=colors[2], label=r"$C_{AB} \propto t^{0.4}$")
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r'$t \ / \ \mathrm{ns}$', fontsize=fsize)
plt.ylabel(r'$C_{AB}$', fontsize=fsize)
plt.legend(fontsize=lsize, frameon=False)


plt.show()