import numpy as np
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

miti = const.eps0.value
angfreqs = np.array(np.logspace(-7, 4, 800))
freqs = angfreqs/(2*np.pi)

base = LayeredSystem(json=r"./Json Models/Enceladus.json")
basecons,aa,bb,cc,dd = base.iterate(angfreqs)

default = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
higher = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
lower = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
thicker = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
thinner = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
more = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")
less = LayeredSystem(json=r"./Json Models/EnceladusWaterLayer.json")

defcon,b,c,d,e = default.iterate(angfreqs)

higher._structure[4][1] = 7500
higher._structure[2][1] = 17500

lower._structure[4][1] = 17500
lower._structure[2][1] = 7500

thicker._structure[2][1] = 10000
thicker._structure[3][1] = 6000
thicker._structure[4][1] = 10000

thinner._structure[2][1] = 12750
thinner._structure[3][1] = 500
thinner._structure[4][1] = 12750

more._structure[3][2] = 2

less._structure[3][2] = 0.5

hc,b,c,d,e = higher.iterate(angfreqs)
lc,b,c,d,e = lower.iterate(angfreqs)
tkc,b,c,d,e = thicker.iterate(angfreqs)
tnc,b,c,d,e = thinner.iterate(angfreqs)
moc,b,c,d,e = more.iterate(angfreqs)
lec,b,c,d,e = less.iterate(angfreqs)

plt.subplot(2,2,1)
plt.grid(True, which="both", ls="--")
plt.plot(freqs, defcon, label = 'Default Water Layer', zorder=8, linewidth=2)
plt.plot(freqs, basecons, color='red', label = 'No Water Layer')
plt.xlim(10**(-3), 500)
plt.ylim(10**(-4.3), 0.4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Water Layer vs No Water Layer')
plt.text(0.02, 0.95, 'a)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.legend()

plt.subplot(2,2,2)
plt.grid(True, which="both", ls="--")
plt.plot(freqs, defcon, label = 'Default Water Layer \n (12.5km deep)', linewidth=2)
plt.plot(freqs, hc, label = 'Closer to Surface \n (7.5 km deep)', linestyle='--', alpha=0.8)
plt.plot(freqs, lc, label='Deeper \n (17.5 km deep)', linestyle='--', alpha=0.8)
plt.plot(freqs, basecons, color='red', label = 'No Water Layer')
plt.xlim(10**(-3), 500)
plt.ylim(10**(-4.3), 0.4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Water Layer Depth')
plt.text(0.02, 0.95, 'b)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.legend()


plt.subplot(2,2,3)
plt.grid(True, which="both", ls="--")
plt.plot(freqs, defcon, label = 'Default Water Layer \n (1 km thick)', zorder=8, linewidth=2)
plt.plot(freqs, tkc, label='Thicker Layer \n (6 km)', linestyle='--', alpha=0.8)
plt.plot(freqs, tnc, label='Thinner Layer \n (500 m)', linestyle='--', alpha=0.8)
plt.plot(freqs, basecons, color='red', label = 'No Water Layer')
plt.xlim(10**(-3), 500)
plt.ylim(10**(-4.3), 0.4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Water Layer Thickness')
plt.text(0.02, 0.95, 'c)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.legend()

plt.subplot(2,2,4)
plt.grid(True, which="both", ls="--")
plt.plot(freqs, defcon, label = 'Default Water Layer \n (1 S/m)', zorder=8, linewidth=2)
plt.plot(freqs, moc, label='More Conductive \n (2 S/m)', zorder=3, alpha=0.8, linestyle='--')
plt.plot(freqs, lec, label='Less Conductive \n (0.5 S/m)', zorder=2, alpha=0.8, linestyle='--')
plt.plot(freqs, basecons, color='red', label = 'No Water Layer')
plt.xlim(10**(-3), 500)
plt.ylim(10**(-4.3), 0.4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]', fontsize=12)
plt.ylabel(r'$\sigma_a$ [S/m]', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.text(0.02, 0.95, 'd)', transform=plt.gca().transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
plt.title('Water Layer Conductivity')
plt.legend()
plt.show()