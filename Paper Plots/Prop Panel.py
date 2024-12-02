import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)

base = LayeredSystem(json=r"./Json Models/Enceladus.json")
basecons,aa,bb,cc,dd = base.iterate(angfreqs)

oconpar = [0.1,0.5,1,2,5]
iconpar = [0.00001,0.0001,0.0005,0.001,0.0015]
ictpar = [10,18,26,35,41]
watpar = [10,18,26,35,41]

oconres = []
for n in oconpar:
    base._structure[1][2] = n
    cons,b,c,d,e = base.iterate(angfreqs)
    oconres.append(cons)

base._structure[1][2] = 1

iconres = []
for n in iconpar:
    base._structure[2][2] = n
    cons,b,c,d,e = base.iterate(angfreqs)
    iconres.append(cons)

base._structure[2][2] = 0.0001

ictres = []
for n in ictpar:
    base._structure[2][1] = n*1000
    cons,b,c,d,e = base.iterate(angfreqs)
    ictres.append(cons)

base._structure[2][1] = 26*1000

watres = []
for n in watpar:
    base._structure[2][1] = n*1000
    base._structure[1][1] = (61-n)*1000
    cons,b,c,d,e = base.iterate(angfreqs)
    watres.append(cons)


plt.subplot(2,2,1)
plt.plot(freqs, oconres[0], label = oconpar[0])
plt.plot(freqs, oconres[1], label = oconpar[1])
plt.plot(freqs, oconres[2], label = oconpar[2])
plt.plot(freqs, oconres[3], label = oconpar[3])
plt.plot(freqs, oconres[4], label = oconpar[4])
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.title(r'Ocean Conductivity [S/m]')
#plt.gca().tick_params(top=True, labeltop=True, right=True, labelright=True)
plt.grid(True, which="both", ls="--")
plt.text(0.02, 0.95, 'a)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.legend()

plt.subplot(2,2,2)
plt.plot(freqs, iconres[0], label=iconpar[0]*1000000)
plt.plot(freqs, iconres[1], label=iconpar[1]*1000000)
plt.plot(freqs, iconres[2], label=iconpar[2]*1000000)
plt.plot(freqs, iconres[3], label=iconpar[3]*1000000)
plt.plot(freqs, iconres[4], label=iconpar[4]*1000000)
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-5.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.title(r'Ice-Shell Conductivity [$\mu$S/m]')
#plt.gca().tick_params(top=True, labeltop=True, right=True, labelright=True)
plt.text(0.02, 0.95, 'b)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.grid(True, which="both", ls="--")
plt.legend()

plt.subplot(2,2,3)
plt.plot(freqs, ictres[0],label =  ictpar[0])
plt.plot(freqs, ictres[1],label = ictpar[1])
plt.plot(freqs, ictres[2],label = ictpar[2])
plt.plot(freqs, ictres[3], label =ictpar[3])
plt.plot(freqs, ictres[4],label = ictpar[4])
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.title('Ice-Shell Thickness [km] \n Constant Ocean Thickness')
#plt.gca().tick_params(top=True, labeltop=True, right=True, labelright=True)
plt.text(0.02, 0.95, 'c)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.grid(True, which="both", ls="--")
plt.legend()

plt.subplot(2,2,4)
plt.plot(freqs,watres[0], label = watpar[0])
plt.plot(freqs,watres[1], label = watpar[1])
plt.plot(freqs,watres[2], label = watpar[2])
plt.plot(freqs,watres[3], label = watpar[3])
plt.plot(freqs,watres[4], label = watpar[4])
plt.xlim(10**(-3), 1000)
plt.ylim(10**(-4.5), 1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.title('Ice-Shell Thickness [km] \n  Constant Total Ocean + Ice Thickness')
#plt.gca().tick_params(top=True, labeltop=True, right=True, labelright=True)
plt.text(0.02, 0.95, 'd)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.grid(True, which="both", ls="--")
plt.legend()

plt.show()