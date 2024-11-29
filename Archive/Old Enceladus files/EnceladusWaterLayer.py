import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

Wl = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
base = LayeredSystem(json=r"./Json Models/Enceladus.json")

miti = const.eps0.value
angfreqs = np.array(np.logspace(-6, 6, 800))
freqs = angfreqs/(2*np.pi)

baseCon, q, prev, int, tanh = base.iterate(angfreqs)
WLconds, q, prev, int, tanh = Wl.iterate(angfreqs)

layerthickness = [0.1,1,5,10,15,20]
layerconds = []
for n in layerthickness:
    Wl._structure[3][2] = n
#    Wl._structure[3][1] = n* 1000
#    Wl._structure[4][1] = (10.5 -(0.5*n)) * 1000
#    Wl._structure[2][1] = (15.5 - (0.5*n)) * 1000
    a,b,c,d,e = Wl.iterate(angfreqs)
    layerconds.append(a)

plt.plot(freqs, baseCon, label='No water layer')
plt.plot(freqs, layerconds[0], label=layerthickness[0])
plt.plot(freqs, layerconds[1], label=layerthickness[1])
plt.plot(freqs, layerconds[2], label=layerthickness[2])
plt.plot(freqs, layerconds[3], label=layerthickness[3])
plt.plot(freqs, layerconds[4], label=layerthickness[4])
plt.xlim(10**(-3), 10000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.title('Water Layer Conductivity [S/m]')
plt.legend()
plt.show()