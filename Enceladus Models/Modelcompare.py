import numpy as np
import matplotlib.pyplot as plt

VHm = np.load(r".\Enceladus Models\VHmost.npy")
VHl = np.load(r".\Enceladus Models\VHleast.npy")
HEMm = np.load(r".\Enceladus Models\HemMost.npy")
HEMl = np.load(r".\Enceladus Models\HemLeast.npy")
Bm = np.load(r".\Enceladus Models\BeutheMost.npy")
Bl = np.load(r".\Enceladus Models\BeutheLeast.npy")

plt.plot(VHm[0], VHm[1], label='Van Hoolst most', linestyle='--', color='blue')
plt.plot(VHl[0], VHl[1], label='Van Hoolst least', linestyle='-.', color='blue')
plt.plot(HEMm[0], HEMm[1], label='Hemmingway most', linestyle='--', color='orange')
plt.plot(HEMl[0], HEMl[1], label='Hemmingway least', linestyle='-.', color='orange')
plt.plot(Bm[0], Bm[1], label='Beuthe most', linestyle='--', color='green')
plt.plot(Bl[0], Bl[1], label='Beuthe least', linestyle='-.', color='green')
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.fill_between(VHm[0], VHm[1], VHl[1], color='blue', alpha=0.3)
plt.fill_between(HEMm[0], HEMm[1], HEMl[1], color='orange', alpha=0.3)
plt.fill_between(Bm[0], Bm[1], Bl[1], color='green', alpha=0.3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)

plt.legend()
plt.show()