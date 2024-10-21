import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from Extractables.ElectromagProps import wavenumber
from Extractables.ElectromagProps import intrin_imp
from Extractables.ElectromagProps import appar_imp
from Extractables.ElectromagProps import appar_resis
from tqdm import tqdm 

miti = const.eps0.value
meab = const.mu0.value


angfreqs = []
intrin = []
squintrin = []
resistivity =[]
params=[150,1000,1560,2500,5000]
list1 = []
list2= []
list3=[]
list4 =[]
list5=[]



for n in params:
    angfreqs = []
    intrin = []
    squintrin = []
    resistivity =[]
    start = 1
    end = 2000000000000000000000
    stepsize = 1

    thickness = n*1000
    conductivity = 10**(-4)
    dieperm = 3.5*miti

    while start< end:
            
            angfreq=start/(10**7)
            angfreqs.append(angfreq)

            intrinsic = intrin_imp(angfreq, meab, conductivity, dieperm)
            intrin.append(np.linalg.norm(intrinsic))
            squareintrin = np.linalg.norm(intrinsic)**2
            squintrin.append(squareintrin)
            resis = appar_resis(intrinsic, angfreq, meab)
            resistivity.append(resis)
            start += stepsize
            stepsize *= 1.05

    if n == params[0]:
        list1.append(resistivity)
    if len(params)>1 and n == params[1]:
        list2.append(resistivity)
    if len(params)>2 and n == params[2]:
        list3.append(resistivity)
    if len(params)>3 and n == params[3]:
        list4.append(resistivity)
    if len(params)>4 and n== params[4]:
        list5.append(resistivity)

    print(angfreqs)

freqs =[]

for entry in angfreqs:
    freqs.append(entry/(2*np.pi))
print(angfreqs)
print(freqs)

condu1=[]
condu2=[]
condu3=[]
condu4=[]
condu5=[]

for entry in list1[0]:
    condu1.append(1/entry)
data = [freqs,condu1]

if len(params)>1:
    for entry in list2[0]:
        condu2.append(1/entry)
if len(params)>2:
    for entry in list3[0]:
        condu3.append(1/entry)
if len(params)>3:
    for entry in list4[0]:
        condu4.append(1/entry)
if len(params)>4:
    for entry in list5[0]:
        condu5.append(1/entry)



plt.plot(freqs, condu1, label=params[0])
if len(params)>1:
    plt.plot(freqs, condu2, label=params[1])
if len(params)>2:
    plt.plot(freqs, condu3, label=params[2])
if len(params)>3:
    plt.plot(freqs, condu4, label=params[3])
if len(params)>4:
    plt.plot(freqs, condu5, label=params[4])
plt.title("Single Ice Layer, Conductivity Varied")
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Resistivity [S/m]')
plt.legend()
plt.show()

