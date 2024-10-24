import numpy as np
import matplotlib.pyplot as plt

OTnocor = np.load(r".\Europa Induction Data\OTnlNoRcorrection.npy")
OT = np.load(r".\Europa Induction Data\OCEANTHICKNESSnolayer.npy")
oceanconductivity = np.load(r".\Europa Induction Data\OCEANCONDUCTIVITYnolayer.npy")
CoreThick = np.load(r".\Europa Induction Data\CoreThickness.npy")
shellthickness = np.load(r".\Europa Induction Data\SHELLTHICKNESSnolayer.npy")
shellthicknocor = np.load(r".\Europa Induction Data\ISTnlNoRcorrection.npy")
ShellPerm = np.load(r".\Europa Induction Data\ISrelPerm.npy")
OceanPerm = np.load(r".\Europa Induction Data\OceanrelPerm.npy")

plt.subplot(2,4,8)
plt.plot(CoreThick[0], CoreThick[1], label='10')
plt.plot(CoreThick[0], CoreThick[2], label='100')
plt.plot(CoreThick[0], CoreThick[3], label='200')
plt.plot(CoreThick[0], CoreThick[4], label='450')
plt.plot(CoreThick[0], CoreThick[5], label='800')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Core Thickness [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-2), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()


plt.subplot(2,4,7)
plt.plot(ShellPerm[0], ShellPerm[1], label='1')
plt.plot(ShellPerm[0], ShellPerm[2], label='3.5')
plt.plot(ShellPerm[0], ShellPerm[3], label='10')
plt.plot(ShellPerm[0], ShellPerm[4], label='20')
plt.plot(ShellPerm[0], ShellPerm[5], label='50')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ice-Shell Relative \n Electric Permittivity")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()


plt.subplot(2,4,6)
plt.plot(OceanPerm[0], OceanPerm[1], label='10')
plt.plot(OceanPerm[0], OceanPerm[2], label='50')
plt.plot(OceanPerm[0], OceanPerm[3], label='85')
plt.plot(OceanPerm[0], OceanPerm[4], label='100')
plt.plot(OceanPerm[0], OceanPerm[5], label='150')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ocean Relative \n Electric Permittivity")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()

plt.subplot(2,4,5)
plt.plot(oceanconductivity[0], oceanconductivity[1], label='0.2')
plt.plot(oceanconductivity[0], oceanconductivity[2], label='0.5')
plt.plot(oceanconductivity[0], oceanconductivity[3], label='1')
plt.plot(oceanconductivity[0], oceanconductivity[4], label='2')
plt.plot(oceanconductivity[0], oceanconductivity[5], label='5')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ocean Conductivity [S/m]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()


plt.subplot(2,4,4)
plt.plot(shellthickness[0], shellthickness[1], label='2')
plt.plot(shellthickness[0], shellthickness[2], label='5')
plt.plot(shellthickness[0], shellthickness[3], label='10')
plt.plot(shellthickness[0], shellthickness[4], label='20')
plt.plot(shellthickness[0], shellthickness[5], label='50')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Shell Thickness \n Radius Corrected [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()


plt.subplot(2,4,3)
plt.plot(shellthicknocor[0], shellthicknocor[1], label='2')
plt.plot(shellthicknocor[0], shellthicknocor[2], label='5')
plt.plot(shellthicknocor[0], shellthicknocor[3], label='10')
plt.plot(shellthicknocor[0], shellthicknocor[4], label='20')
plt.plot(shellthicknocor[0], shellthicknocor[5], label='50')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ice-Shell Thickness\n Radius NOT Corrected [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-4), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()

plt.subplot(2,4,2)
plt.plot(OT[0], OT[1], label='40')
plt.plot(OT[0], OT[2], label='60')
plt.plot(OT[0], OT[3], label='100')
plt.plot(OT[0], OT[4], label='160')
plt.plot(OT[0], OT[5], label='250')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ocean Thickness\n Radius Corrected [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-2), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()


plt.subplot(2,4,1)
plt.plot(OTnocor[0], OTnocor[1], label='40')
plt.plot(OTnocor[0], OTnocor[2], label='60')
plt.plot(OTnocor[0], OTnocor[3], label='100')
plt.plot(OTnocor[0], OTnocor[4], label='160')
plt.plot(OTnocor[0], OTnocor[5], label='250')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Apparent Conductivity [S/m]')
plt.title("Ocean Thickness\n Radius NOT Corrected [km]")
plt.xlim(10**(-6), 10)
plt.ylim(10**(-2), 10)
plt.axvspan(10**(-2.5), 100, color='lightgrey')
plt.axvline(x=4.3*10**(-7), color='black', linestyle='--')
plt.axvline(x=3.3*10**(-6), color='black', linestyle='--')
plt.axvline(x=2.5*10**(-5), color='black', linestyle='--')
plt.axvline(x=4.9*10**(-5), color='black', linestyle='--')
plt.axvline(x=7.4*10**(-5), color='black', linestyle='--')
plt.legend()

plt.suptitle("Europa Induction Curves with No Water Layer", fontsize= 22)
plt.show()