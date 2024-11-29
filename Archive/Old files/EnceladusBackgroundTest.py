import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)
model = LayeredSystem(json=r"./Json Models/EnceladusBackground.json")

params = [0.0001,0.001,0.01,0.1,1]
condu = []

for n in params:
    model._structure[0][2] = n
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    condu.append(conds)

plt.plot(freqs, condu[0], label=params[0])
plt.plot(freqs, condu[1], label=params[1])
plt.plot(freqs, condu[2], label=params[2])
plt.plot(freqs, condu[3], label=params[3])
plt.plot(freqs, condu[4], label=params[4])
plt.yscale('log')
plt.xscale('log')
plt.xlim(10**(-6), 1000)
plt.ylim(10**(-5), 10)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.legend()
plt.title('Background Ocean Conductivtiy [S/m]')
plt.show()
