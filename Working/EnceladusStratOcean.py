import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 5, 800))
freqs = angfreqs/(2*np.pi)
base = LayeredSystem(json=r"./Json Models/Enceladus.json")
test = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")

baseCon, q, prev, int, tanh = base.iterate(angfreqs)
testconds, tq, tprev, tint, ttanh = test.iterate(angfreqs)

plt.plot(freqs, baseCon, label='Not Stratified')
plt.plot(freqs, testconds, label='Stratified')
plt.xlim(10**(-3), 10000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title('Water Layer Conductivity [S/m]')
plt.legend()
plt.show()
