import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt

meab = const.mu0.value
cond = 1
omega = []
depths = []

start = 1
end = 2000000000000
stepsize = 1

while start < end:

    freq = start/ (10**7)
    depth = np.sqrt(2/(meab*cond*freq))
    depths.append(depth/1000)

    omega.append(freq)

    start += stepsize
    stepsize *= 1.05

freqs = []

for entry in omega:
    freqs.append(entry/(2*np.pi))

plt.plot(freqs, depths)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Skin Depth [km]')
plt.xscale('log')
plt.yscale('log')
plt.show()