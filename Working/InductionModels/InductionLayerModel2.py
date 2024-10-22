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
params = [10]
list1 = []
list2= []
list3=[]
list4 =[]
list5=[]

for n in params:
    resistivities=[]
    angfreqs=[]
    start = 1
    end = 2000000000
    stepsize = 1
    while start < end:
    
        angfreq=start/(10**7)
        angfreqs.append(angfreq)

        magperm1 = 5000*meab                  #Innermost layer     iron core
        thick1 = 450000
        cond1 = 10**6
        dieperm1 = 10000*miti

        wave1 = wavenumber(angfreq, magperm1, cond1, dieperm1)
        intrin1 = intrin_imp(angfreq, magperm1, cond1, dieperm1)
        apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)
        apparimps.append(apparimp1)

        magperm2 = 1*meab                 #Second innermost layer    silicate rock
        thick2 = 1000*1000
        cond2 = 10**(-3)
        dieperm2 = 5*miti

        wave2 = wavenumber(angfreq, magperm2, cond2, dieperm2)
        intrin2 = intrin_imp(angfreq, magperm2, cond2, dieperm2)
        apparimp2 = appar_imp(wave2, intrin2, apparimp1, thick2)
        apparimps.append(apparimp2)


        magperm3 = 1*meab                  #Third layer                  water
        thick3 = 100*1000
        cond3 = 1
        dieperm3 =85*miti

        wave3 = wavenumber(angfreq, magperm3, cond3, dieperm3)
        intrin3 = intrin_imp(angfreq, magperm3, cond3, dieperm3)
        apparimp3 = appar_imp(wave3, intrin3, apparimp2, thick3)
        apparimps.append(apparimp3)

        magperm4 = 1*meab                  #Fourth layer                     ice
        thick4 = n*1000
        cond4 = 10**(-4)
        dieperm4 = 3.5*miti

        wave4 = wavenumber(angfreq, magperm4, cond4, dieperm4)
        intrin4 = intrin_imp(angfreq, magperm4, cond4, dieperm4)
        apparimp4 = appar_imp(wave4, intrin4, apparimp3, thick4)
        apparimps.append(apparimp4)

        """

        magperm5 = 1*meab                   #Fifth layer                      water
        thick5 = 2000
        cond5 = 1
        dieperm5 = 85*miti

        wave5 = wavenumber(angfreq, magperm5, cond5, dieperm5)
        intrin5 = intrin_imp(angfreq, magperm5, cond5, dieperm5)
        apparimp5 = appar_imp(wave5, intrin5, apparimp4, thick5)
        apparimps.append(apparimp5)

        magperm6 = 1*meab                  #Sixth layer                        ice
        thick6 = 5000
        cond6 = 10**(-4)
        dieperm6 = 3.5*miti

        wave6 = wavenumber(angfreq, magperm6, cond6, dieperm6)
        intrin6 = intrin_imp(angfreq, magperm6, cond6, dieperm6)
        apparimp6 = appar_imp(wave6, intrin6, apparimp5, thick6)
        apparimps.append(apparimp6)
        """

        resistivity = appar_resis(apparimps[-1], angfreq, magperm4)
        resistivities.append(resistivity)


        start += stepsize
        stepsize *= 1.05
    
    if n == params[0]:
        list1.append(resistivities)
    if len(params)>1 and n == params[1]:
        list2.append(resistivities)
    if len(params)>2 and n == params[2]:
        list3.append(resistivities)
    if len(params)>3 and n == params[3]:
        list4.append(resistivities)
    if len(params)>4 and n== params[4]:
        list5.append(resistivities)

    

condu1 = []
condu2 = []
condu3 = []
condu4 = []
condu5 = []

freqs =[]

for entry in angfreqs:
    freqs.append(entry/(2*np.pi))


for entry in list1[0]:
    condu1.append(1/entry)
data = [freqs,condu1]

if len(params)>1:
    for entry in list2[0]:
        condu2.append(1/entry)
    data = [freqs,condu1,condu2]
if len(params)>2:
    for entry in list3[0]:
        condu3.append(1/entry)
    data = [freqs, condu1,condu2,condu3]
if len(params)>3:
    for entry in list4[0]:
        condu4.append(1/entry)
    data = [freqs,condu1,condu2,condu3,condu4]
if len(params)>4:
    for entry in list5[0]:
        condu5.append(1/entry)
    data = [freqs,condu1,condu2,condu3,condu4,condu5]


plt.plot(freqs, condu1, label=params[0])
if len(params)>1:
    plt.plot(freqs, condu2, label=params[1])
if len(params)>2:
    plt.plot(freqs, condu3, label=params[2])
if len(params)>3:
    plt.plot(freqs, condu4, label=params[3])
if len(params)>4:
    plt.plot(freqs, condu5, label=params[4])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Core Thickness [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()

plt.show()


np.save(r".\Europa Models\BaseCurve.npy", data)