import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 5, 800))
freqs = angfreqs/(2*np.pi)
base = LayeredSystem(json=r"./Json Models/Enceladus.json")
test1 = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")
test2 = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")
test = LayeredSystem(json=r"./Json Models/EnceladusOceanStrat.json")

baseCon, q, prev, int, tanh = base.iterate(angfreqs)
#testconds, tq, tprev, tint, ttanh = test.iterate(angfreqs)

iterconds = []
params = [0.01,0.1,1,5,20]

#for n in params:
#    test._structure[1][1] = n*1000
#    test._structure[2][1] = (35-n)*1000
#    test._structure[1][2] = n
#    testconds, tq, tprev, tint, ttanh  = test.iterate(freqs)
#    iterconds.append(testconds)

test1._structure[1][1] = 5000
test1._structure[1][2] = 20
test1._structure[2][1] = 30000
test1._structure[2][2] = 0.0001

test2._structure[1][1] = 25000
test2._structure[1][2] = 1
test2._structure[2][1] = 10000
test2._structure[2][2] = 0.01

test1cond, b,c,d,e = test1.iterate(freqs)
test2cond, f,g,h,i = test2.iterate(freqs)

#plt.plot(freqs, baseCon, label='Not Stratified')
#plt.plot(freqs, iterconds[0], label=params[0])
#plt.plot(freqs, iterconds[1], label=params[1])
#plt.plot(freqs, iterconds[2], label=params[2])
#plt.plot(freqs, iterconds[3], label=params[3])
#plt.plot(freqs, iterconds[4], label=params[4])

plt.plot(freqs, test1cond, label='Bigger Inverse Layer')
plt.plot(freqs, test2cond, label='Smaller Inverese Layer')
plt.xlim(10**(-3), 10000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title('Situational Comparison')
plt.legend()
plt.show()
