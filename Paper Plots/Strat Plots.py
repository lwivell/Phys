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

stratper =100 - (stratcon/basecon)*100
gradper =100 - (gradcon/basecon)*100

plt.subplot(3,1,1)
plt.grid()
plt.plot(freqs, stratcon, label='Stratified Ocean')
plt.plot(freqs, gradcon, label='Ocean Conductivity Gradient')
plt.plot(freqs, basecon, label='Constant Ocean Conductivity')
plt.xlim(10**(-3), 100)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.text(0.02, 0.95, 'a)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.xticks(fontsize=0)
plt.yticks(fontsize=10)
plt.legend()

plt.subplot(3,1,2)
plt.grid()
plt.plot(freqs, stratdiff, label='Stratified Ocean')
plt.plot(freqs, graddiff, label='Ocean Conductivity Gradient')
plt.axhline(y=0, color='gray', linestyle='--', label = 'Constant Ocean Conductivity')
plt.xlim(10**(-3), 100)
plt.ylim(-0.019, 0.002)
plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$\Delta$ $\sigma_a$ [S/m]', fontsize=12)
plt.text(0.02, 0.85, 'b)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.xticks(fontsize=0)
plt.yticks(fontsize=10)
plt.legend()


plt.subplot(3,1,3)
plt.grid()
plt.plot(freqs, stratper, label='Stratified Ocean')
plt.plot(freqs, gradper, label='Ocean Conductivity Gradient')
#plt.plot(freqs, basecon)
plt.xlim(10**(-3), 100)
#plt.ylim(10**(-4.5), 1)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel('Percentage Difference [%]', fontsize=12)
plt.axhline(y=0, color='gray', linestyle='--', label = 'Constant Ocean Conductivity')
plt.text(0.02, 0.95, 'c)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
#plt.title('Situational Comparison')
plt.legend()
plt.show()
