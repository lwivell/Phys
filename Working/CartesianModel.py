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

apparimps = []
data = []
resistivities=[]
angfreqs=[]
for jumps in tqdm(np.arange(1, 10**7, 100000)):

    angfreq=jumps/(10**6)
    angfreqs.append(angfreq)

    magperm1 = 5000*meab                  #Innermost layer     iron core
    thick1 = 450000
    cond1 = 10**6
    dieperm1 = 1*miti

    wave1 = wavenumber(angfreq, magperm1, cond1, dieperm1)
    intrin1 = intrin_imp(angfreq, magperm1, cond1, dieperm1)
    apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)
    apparimps.append(apparimp1)

    magperm2 = 1*meab                 #Second innermost layer    silicate rock
    thick2 = 1000000
    cond2 = 10**(-3)
    dieperm2 = 1*miti

    wave2 = wavenumber(angfreq, magperm2, cond2, dieperm2)
    intrin2 = intrin_imp(angfreq, magperm2, cond2, dieperm2)
    apparimp2 = appar_imp(wave2, intrin2, apparimp1, thick2)
    apparimps.append(apparimp2)


    magperm3 = 1*meab                  #Third layer                  water
    thick3 = 100000
    cond3 = 1
    dieperm3 = 1*miti

    wave3 = wavenumber(angfreq, magperm3, cond3, dieperm3)
    intrin3 = intrin_imp(angfreq, magperm3, cond3, dieperm3)
    apparimp3 = appar_imp(wave3, intrin3, apparimp2, thick3)
    apparimps.append(apparimp3)

    magperm4 = 1*meab                  #Fourth layer                     ice
    thick4 = 5000
    cond4 = 10**(-4)
    dieperm4 = 1*miti

    wave4 = wavenumber(angfreq, magperm4, cond4, dieperm4)
    intrin4 = intrin_imp(angfreq, magperm4, cond4, dieperm4)
    apparimp4 = appar_imp(wave4, intrin4, apparimp3, thick4)
    apparimps.append(apparimp4)

    magperm5 = 1*meab                   #Fifth layer                      water
    thick5 = 3000
    cond5 = 1
    dieperm5 = 1*miti

    wave5 = wavenumber(angfreq, magperm5, cond5, dieperm5)
    intrin5 = intrin_imp(angfreq, magperm5, cond5, dieperm5)
    apparimp5 = appar_imp(wave5, intrin5, apparimp4, thick5)
    apparimps.append(apparimp5)

    magperm6 = 1*meab                  #Sixth layer                        ice
    thick6 = 2000
    cond6 = 10**(-4)
    dieperm6 = 1*miti

    wave6 = wavenumber(angfreq, magperm6, cond6, dieperm6)
    intrin6 = intrin_imp(angfreq, magperm6, cond6, dieperm6)
    apparimp6 = appar_imp(wave1, intrin6, apparimp5, thick6)
    apparimps.append(apparimp6)

    resistivity = appar_resis(apparimps[-1], angfreq, magperm6)
    resistivities.append(resistivity)
    data.append([angfreq, resistivity])

conductivities=[]
for entry in resistivities:
    conductivities.append(1/entry)

print(conductivities)
print(angfreqs)
plt.subplot(2,1,1)
plt.plot(angfreqs, conductivities)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Angular frequencies [rads]')
plt.ylabel('Apparent Conductivity')

plt.subplot(2,1,2)
plt.plot(angfreqs, resistivities)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Angular frequencies [rads]')
plt.ylabel('Apparent Conductivity')

plt.show()


np.save("resisData.npy", data)


