import numpy as np 
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-9, 7, 800))
freqs = angfreqs/(2*np.pi)

grad = LayeredSystem(json=r"./Json Models/Gradient.json")

conds, b,c,d,e = grad.iterate(angfreqs)


plt.plot(freqs, conds)
#plt.xlim(10**(-3), 10000)
#plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title('Nineteen layers')
plt.legend()
plt.show()