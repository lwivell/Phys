import numpy as np 
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-9, 7, 800))
freqs = angfreqs/(2*np.pi)

grad = LayeredSystem(json=r"./Json Models/EnceladusOceanGrad.json")
strat = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")
base = LayeredSystem(json=r"./Json Models/Enceladus.json")

conds, b,c,d,e = grad.iterate(angfreqs)
stratconds,f,g,h,i = strat.iterate(angfreqs)
baseconds, j,k,l,m = base.iterate(angfreqs)

plt.plot(freqs, conds, label='Ocean Conductivity Gradient')
plt.plot(freqs, stratconds, label='Stratified Ocean')
plt.plot(freqs, baseconds, label='Ocean Conductivity Constant')
plt.xlim(10**(-3), 10)
plt.ylim(10**(-4.5), 0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title('Comparison of Model Regimes')
plt.legend()
plt.show()