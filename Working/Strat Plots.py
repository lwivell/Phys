import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)
strat = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")
base = LayeredSystem(json=r"./Json Models/Enceladus.json")
grad = LayeredSystem(json=r"./Json Models/EnceladusOceanGrad.json")


basecon, q, prev, int, tanh = base.iterate(angfreqs)
stratcon, a,b,c,d = strat.iterate(angfreqs)
gradcon, e,f,g,h, = grad.iterate(angfreqs)

stratdiff = stratcon - basecon
graddiff = gradcon - basecon

plt.plot(freqs, stratdiff, label='Stratified Ocean')
plt.plot(freqs, graddiff, label='Ocean Conductivity Gradient')
plt.axhline(y=0, color='gray', linestyle='--', label = 'No Difference')
plt.xlim(10**(-3), 1)
plt.ylim(-0.02, 0.002)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\Delta$ $\sigma_a$ [S/m]', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend()
plt.show()

plt.plot(freqs, stratcon, label='Stratified')
plt.plot(freqs, gradcon, label='Gradient')
plt.plot(freqs, basecon)
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
#plt.title('Situational Comparison')
plt.legend()
plt.show()
