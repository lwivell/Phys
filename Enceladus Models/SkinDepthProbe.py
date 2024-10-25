import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

data = np.load(r".\Enceladus Models\BaseCurve.npy")
meab = const.mu0.value


freqs = data[0]
cond = data[1]
depths = []

denoms = [a*b for a,b in zip(freqs,cond)]

for n in denoms:
    z= np.sqrt(2/(meab*n))
    depths.append(z)

scale = []
for n in cond:
    t = n*meab*(252000**2)
    scale.append(t)

Periods = []
for n in freqs:
    P = (2*np.pi)/n
    Periods.append(P)

kmdepths = []
for entry in depths:
    kmdepths.append(entry/1000)

zgold = []
for entry in kmdepths:
    zgold.append(np.linalg.norm(entry/252)**2)

ratio = [a/b for a,b in zip(Periods, scale)]

tgold= []
for entry in ratio:
    tgold.append((entry**2)/np.pi)



plt.subplot(1,2,1)
plt.plot(freqs, scale, label='Diffusion Time')
plt.xscale('log')
plt.xlim(10**(-6), 10)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Diffusion Time [s]')

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

plt.plot(freqs, zgold, label='Length Scale Ratio')
plt.plot(freqs, tgold, label='Time Scale Ratio')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
#plt.xlim(10**(-6), 10)
plt.axhline(y=1, color='black', linestyle='--')
plt.legend()
plt.show()

plt.plot(freqs, kmdepths)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Skin Depth [km]')
plt.grid(True, which="both", ls="--")
plt.show()

np.save(r'.\Enceladus Models\apparskin.npy', [freqs, kmdepths])