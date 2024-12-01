import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)

base = LayeredSystem(json=r"./Json Models/Enceladus.json")
basecons,aa,bb,cc,dd = base.iterate(angfreqs)










plt.subplot(2,2,1)
plt.plot(freqs, basecon)
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.title('Situational Comparison')
plt.legend()









plt.show()