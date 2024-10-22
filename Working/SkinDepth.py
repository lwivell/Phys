import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

data = np.load(r".\Europa Models\BaseCurve.npy")
meab = const.mu0.value


freqs = data[0]
cond = data[1]
depths = []

denoms = [a*b for a,b in zip(freqs,cond)]

for n in denoms:
    z= np.sqrt(2/(meab*n))
    depths.append(z)


kmdepths =[]
for entry in depths:
    kmdepths.append(entry/1000)

plt.subplot(1,2,1)
plt.plot(freqs, cond)
plt.yscale('log')
plt.xscale('log')
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')

plt.subplot(1,2,2)
plt.plot(freqs, kmdepths)
plt.axhline(y = 1110, color='orange', linestyle='--', label='Depth with Mantle')
plt.axhline(y=1560, color='red', linestyle='--', label='Depth with Core')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Skin Depth [km]')
plt.xlim(10**(-6), 10)
plt.legend()
plt.show()