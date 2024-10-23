import numpy as np
import matplotlib.pyplot as plt

whole = np.load(r".\Enceladus Models\whole.npy")
tophalf = np.load(r".\Enceladus Models\BaseCurve.npy")
bottomhalf = np.load(r".\Enceladus Models\bottomhalf.npy")

plt.plot(whole[0], whole[1], label = 'Whole Body')
plt.plot(tophalf[0], tophalf[1], label='Top Half of Body', linestyle = '--')
plt.plot(bottomhalf[0], bottomhalf[1], label='Bottom Half of Body', linestyle='--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
#plt.title('Mantle Thickness [km]')
#plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.legend()

plt.show()