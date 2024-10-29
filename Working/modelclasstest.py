import numpy as np
from Extractables.InductionModel import LayeredSystem
import matplotlib.pyplot as plt

test1 = LayeredSystem(json=r"./Json Models/Enceladus.json")
test2 = LayeredSystem(json=r"./Json Models/Europa.json")

angfreqs = np.array(np.logspace(-7,2,400))
conds1 = test1.iterate(angfreqs)
conds2 = test2.iterate(angfreqs)

freqs = angfreqs/(2*np.pi)

plt.plot(freqs, conds1)
plt.plot(freqs, conds2)
plt.xscale('log')
plt.yscale('log')
plt.show()