import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

"""
For iterating, important to understand the class storage structure

model._structure[n][x]

n=0 is the innermost layer, moving outwards

x=0 is the index of the layer
x=1 is the thickness of the layer
x=2 is the conductivity of the layer
x=3 is the electric permittivity of the layer

Use these indexes in the for loop to iterate through different trial values
"""

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)
model = LayeredSystem(json=r"./Json Models/Enceladus.json")

OceanThickness = [20,28,35,42,50]                  # Iterates corrected ocean thickness
OTresponse=[]
for n in OceanThickness:
    model._structure[1][1] = n*1000
    model._structure[2][1] = (61-n)*1000
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    OTresponse.append(conds)

model._structure[1][1] = 35000
model._structure[2][1] = 26000

CoreThickness = [150,175,190,200,225]              #Iterates core mantle thickness
CTresponse = []
for n in CoreThickness:
    model._structure[0][1] = n*1000
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    CTresponse.append(conds)   

model._structure[0][1] = 190000

IceCond = [0.00001, 0.0001, 0.001, 0.01,0.1]      # Iterates ice shell conductivity
ICresponse = []
for n in IceCond:
    model._structure[2][2] = n
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    ICresponse.append(conds)

model._structure[2][2] = 0.0001

OceanCond = [0.01,0.1,1,5,10]                       #Iterates ocean Conductivity
OCresponse = []
for n in OceanCond:
    model._structure[1][2] = n
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    OCresponse.append(conds)

model._structure[1][2] = 1

CoreCond = [0.00001, 0.0001, 0.001, 0.01,0.1]        # Iterates core conductivity
CCresponse = []
for n in CoreCond:
    model._structure[0][2] = n
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    CCresponse.append(conds)    

model._structure[0][2] = 0.001

CorePerm = [1,3.5,5,20,85]                             # Iterates core permittivities
CPresponse = []
for n in CorePerm:
    model._structure[0][3] = n*miti
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    CPresponse.append(conds)  

model._structure[0][3] = 5 * miti

OceanPerm = [1,3.5,5,20,85]                            # Iterates ocean permittivities
OPresponse = []
for n in CorePerm:
    model._structure[1][3] = n*miti
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    OPresponse.append(conds)

model._structure[1][3] = 85 * miti

IcePerm = [1,3.5,5,20,85]                              # Iterates ice permittivities
IPresponse = []
for n in CorePerm:
    model._structure[2][3] = n*miti
    conds, q, prev, int, tanh = model.iterate(angfreqs)
    IPresponse.append(conds)  

model._structure[2][3] = 3.5 * miti

plt.plot(freqs, OCresponse[0], label=OceanCond[0])
plt.plot(freqs, OCresponse[1], label=OceanCond[1])
plt.plot(freqs, OCresponse[2], label=OceanCond[2])
plt.plot(freqs, OCresponse[3], label=OceanCond[3])
plt.plot(freqs, OCresponse[4], label=OceanCond[4])
plt.yscale('log')
plt.xscale('log')
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-5), 10)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.legend()
plt.title('Ocean Conductivity [S/m]')
plt.show()