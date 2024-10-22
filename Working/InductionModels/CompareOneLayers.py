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

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

thickness = 10*1000
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

angfreqs = []
intrin = []
squintrin = []
resistivity =[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

thickness = 100*1000
conductivity = 1
dieperm = 85*miti

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

watercondu =[]
for entry in resistivity:
      watercondu.append(1/entry)

angfreqs = []
intrin = []
squintrin = []
resistivity =[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

thickness = 1000*1000
conductivity = 10**(-3)
dieperm = 5*miti

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

silicondu =[]
for entry in resistivity:
      silicondu.append(1/entry)

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

apparimps=[]
resistivities=[]

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        magperm1 = meab                  #Water
        thick1 = 100*1000
        cond1 = 1
        dieperm1 = 85*miti

        wave1 = wavenumber(angfreq, magperm1, cond1, dieperm1)
        intrin1 = intrin_imp(angfreq, magperm1, cond1, dieperm1)
        apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)
        apparimps.append(apparimp1)

        magperm2 = 1*meab                 #Ice
        thick2 = 10*1000
        cond2 = 10**(-4)
        dieperm2 = 3.5*miti

        wave2 = wavenumber(angfreq, magperm2, cond2, dieperm2)
        intrin2 = intrin_imp(angfreq, magperm2, cond2, dieperm2)
        apparimp2 = appar_imp(wave2, intrin2, apparimp1, thick2)
        apparimps.append(apparimp2)

        start += stepsize
        stepsize *= 1.05
        resistivity = appar_resis(apparimps[-1], angfreq, meab)
        resistivities.append(resistivity)


twocondu1=[]
for entry in resistivities:
    twocondu1.append(1/entry)

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

apparimps=[]
resistivities=[]

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        magperm1 = meab                  #silicate 
        thick1 = 1000*1000
        cond1 = 10**(-3)
        dieperm1 = 5*miti

        wave1 = wavenumber(angfreq, magperm1, cond1, dieperm1)
        intrin1 = intrin_imp(angfreq, magperm1, cond1, dieperm1)
        apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)
        apparimps.append(apparimp1)

        magperm2 = 1*meab                 #water
        thick2 = 100*1000
        cond2 = 1
        dieperm2 = 85*miti

        wave2 = wavenumber(angfreq, magperm2, cond2, dieperm2)
        intrin2 = intrin_imp(angfreq, magperm2, cond2, dieperm2)
        apparimp2 = appar_imp(wave2, intrin2, apparimp1, thick2)
        apparimps.append(apparimp2)

        start += stepsize
        stepsize *= 1.05
        resistivity = appar_resis(apparimps[-1], angfreq, meab)
        resistivities.append(resistivity)


twocondu2=[]
for entry in resistivities:
    twocondu2.append(1/entry)



resistivithree =[]

start = 1
end = 2000000000000000000000000000000000
stepsize = 1

while start< end:
        
        angfreq=start/(10**13)
        angfreqs.append(angfreq)

        magperm1 = meab                  #Innermost layer     Silicate
        thick1 = 1000*1000
        cond1 = 10**(-3)
        dieperm1 = 5*miti

        wave1 = wavenumber(angfreq, magperm1, cond1, dieperm1)
        intrin1 = intrin_imp(angfreq, magperm1, cond1, dieperm1)
        apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)
        apparimps.append(apparimp1)

        magperm2 = meab                 #Second innermost layer    water
        thick2 = 100*1000
        cond2 = 1
        dieperm2 = 85*miti

        wave2 = wavenumber(angfreq, magperm2, cond2, dieperm2)
        intrin2 = intrin_imp(angfreq, magperm2, cond2, dieperm2)
        apparimp2 = appar_imp(wave2, intrin2, apparimp1, thick2)
        apparimps.append(apparimp2)


        magperm3 = meab                  #Third layer             ice
        thick3 = 10*1000
        cond3 = 10**(-4)
        dieperm3 =3.5*miti

        wave3 = wavenumber(angfreq, magperm3, cond3, dieperm3)
        intrin3 = intrin_imp(angfreq, magperm3, cond3, dieperm3)
        apparimp3 = appar_imp(wave3, intrin3, apparimp2, thick3)
        apparimps.append(apparimp3)

        resistivity = appar_resis(apparimps[-1], angfreq, meab)
        resistivithree.append(resistivity)


        start += stepsize
        stepsize *= 1.05



threecondu=[]
for entry in resistivithree:
    threecondu.append(1/entry)


plt.title("Components of an Induction Curve")
plt.plot(freqs, silicondu, label ='Silicate Layer', linestyle='--')
plt.plot(freqs, watercondu, label='Water Layer', linestyle='--')
plt.plot(freqs, icecondu, label='Ice Layer', linestyle='--')
plt.plot(freqs, twocondu1, label='Ice + Water Layers', linestyle='-.')
plt.plot(freqs, twocondu2, label='Silicate + Water', linestyle='-.')
plt.plot(freqs, threecondu, label='Complete Model')
plt.yscale('log')
plt.xscale('log')
plt.xlim(10**(-8), 1000)
plt.ylim(10**(-5), 10)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.legend()
plt.show()

