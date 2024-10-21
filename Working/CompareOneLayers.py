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

list1 = []
list2= []
list3=[]
list4 =[]
list5=[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

thickness = 1560*1000
conductivity = 10**(-4)
dieperm = 3.5*miti

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        intrinsic = intrin_imp(angfreq, meab, conductivity, dieperm)
        intrin.append(np.linalg.norm(intrinsic))
        squareintrin = np.linalg.norm(intrinsic)**2
        squintrin.append(squareintrin)
        resis = appar_resis(intrinsic, angfreq, meab)
        resistivity.append(resis)
        start += stepsize
        stepsize *= 1.05

freqs =[]
for entry in angfreqs:
    freqs.append(entry/(2*np.pi))
icecondu=[]
for entry in resistivity:
    icecondu.append(1/entry)

conductivity = 10**(-3)
dieperm = 5*miti

resistivity =[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        intrinsic = intrin_imp(angfreq, meab, conductivity, dieperm)
        intrin.append(np.linalg.norm(intrinsic))
        squareintrin = np.linalg.norm(intrinsic)**2
        squintrin.append(squareintrin)
        resis = appar_resis(intrinsic, angfreq, meab)
        resistivity.append(resis)
        start += stepsize
        stepsize *= 1.05



silicondu=[]
for entry in resistivity:
    silicondu.append(1/entry)

conductivity = 1
dieperm = 85*miti

resistivity =[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        intrinsic = intrin_imp(angfreq, meab, conductivity, dieperm)
        intrin.append(np.linalg.norm(intrinsic))
        squareintrin = np.linalg.norm(intrinsic)**2
        squintrin.append(squareintrin)
        resis = appar_resis(intrinsic, angfreq, meab)
        resistivity.append(resis)
        start += stepsize
        stepsize *= 1.05



watercondu=[]
for entry in resistivity:
    watercondu.append(1/entry)

plt.title("Single Layer Models")
plt.plot(freqs, watercondu, label='Water')
plt.plot(freqs, icecondu, label='Ice')
plt.plot(freqs, silicondu, label='Silicate')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Conductivity [S/m]')
plt.legend()
plt.show()

