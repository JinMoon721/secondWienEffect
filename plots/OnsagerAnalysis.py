# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 12:48:02 2026

@author: whisk
"""

from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE

from scipy.special import i1
qqrd2e =332.06371 ## e^2/4 pi epsilon0 in kcal/mol A unit
kBT = 0.592 ## kcal/mol
dielectric=69.565
#dielectric=19.202

Bjerrum = qqrd2e/kBT/dielectric ## A unit

field = 50 ## mV/A
kcal2mev = 25.7/kBT
qE = field/kcal2mev ## kcal/molA 

elength=kBT/qE ## A unit

x=Bjerrum/elength/2 ## factor of 2 from Onsager's definition


typicalLength = np.sqrt( 20 * 0.18) * 10 
print(typicalLength)

prev=np.exp( 4.34/kcal2mev * 18.9 /kBT)

print(np.log(prev))


## what if this is barrier crossing event
length=(9.1-5.3)/4
length=0.2
deltaF = field/kcal2mev * length / kBT
print(deltaF*kBT)
print("1D Barrier Crossing : %f" % (np.exp(deltaF)))


## heat bound
D=(5.9+7) ## nm^2/ns
tau=12 ## ps

D=(0.6+0.9)

drift= D*100 * 50 /25 ## A/ns
dist=drift * tau /1000 ## A
print("drift velocity= %lf" % (drift) ) 
print("dist = %lf" % (dist))

heatbound = 1/kBT**2 * (field/kcal2mev)**2 * (D*100 /1000) * tau /4
print("heatbound : %f" % (np.exp(heatbound)))

ratio = i1(np.sqrt(8*x)) / np.sqrt(2*x)
print("F(x) = %f" % (ratio) )
conc=0.5
#k0= 1.753 * 10**(-5)
bq=0.401
k0=bq**2 * conc / (1-bq)
print("k0=%lf" % (k0))

k=k0 *ratio
print(np.sqrt(ratio))

newratio = (np.sqrt( k**2 + 4*k*conc) - k) / (np.sqrt(k0**2 + 4*k0*conc) - k0)
print("non-strong case: %lf" % (newratio))


field = np.linspace(0, 50, num=100)
#dielectric=20
qE = field/kcal2mev ## kcal/molA 
elength=kBT/qE ## A unit
Bjerrum = qqrd2e/kBT/dielectric ## A unit
x=Bjerrum/elength/2 ## factor of 2 from Onsager's definition


ratio = i1(np.sqrt(8*x)) / np.sqrt(2*x)
k=k0*ratio
newratio = (np.sqrt( k**2 + 4*k*conc) - k) / (np.sqrt(k0**2 + 4*k0*conc) - k0)

david = 1/2* ( np.sqrt(ratio*(4/k0+ratio)) - ratio)
fig, ax = plt.subplots(figsize=(5,5))
plt.plot(field, newratio)
plt.plot(field, field*0 + 2.4018)
#plt.plot(field, david)
#plt.plot(field, np.exp(field/60)/1.13, '--', lw=3)

#plt.yscale("log")
plt.show()