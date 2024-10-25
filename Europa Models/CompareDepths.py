import matplotlib.pyplot as plt
import numpy as np

effect = np.load(r".\Europa Models\EffectiveSkinDepth.npy")
appar = np.load(r".\Europa Models\ApparentSkinDepth.npy")

plt.plot(effect[0], effect[1], label = 'Effective Skin Depth')
plt.plot(appar[0], appar[1], label='Apparent Skin Depth')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequncy [Hz]')
plt.ylabel('Skin Depth [km]')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()