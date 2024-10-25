import numpy as np
import matplotlib.pyplot as plt

appar =np.load(r'.\Enceladus Models\apparskin.npy')
effect = np.load(r'.\Enceladus Models\EffectiveSkinDepth.npy')

plt.plot(effect[0], effect[1], label = 'Effective Skin Depth')
plt.plot(appar[0], appar[1], label='Apparent Skin Depth')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequncy [Hz]')
plt.ylabel('Skin Depth [km]')
plt.axhline(y=251, color='black', linestyle='--')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()