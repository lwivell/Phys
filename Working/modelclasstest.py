import numpy as np
from Extractables.InductionModel import LayeredSystem
import matplotlib.pyplot as plt

test1 = LayeredSystem(json=r"./Json Models/Enceladus.json")
test2 = LayeredSystem(json=r"./Json Models/Europa.json")

angfreqs = np.array(np.logspace(-7,2,400))
conds1, q1 = test1.iterate(angfreqs)
conds2, q2 = test2.iterate(angfreqs)

freqs = angfreqs/(2*np.pi)

plt.plot(freqs, conds1)
plt.plot(freqs, conds2)
plt.xscale('log')
plt.yscale('log')
plt.show()

fig, ax1 = plt.subplots()

line1 = ax1.plot(freqs, conds1, label = 'Apparent Conductivity', color='orange')
ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Apparent Conductivity [S/m]')
ax1.set_yscale('log')
ax1.set_xscale('log')


ax2 = ax1.twinx()


line2 = ax2.plot(freqs, q1, label='Q-function', color='blue')
ax2.set_ylabel('Q-function value')
ax2.set_yscale('log')

lines = line1 + line2
labels = [ line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper left')
plt.show()