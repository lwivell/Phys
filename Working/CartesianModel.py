import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from Extractables.ElectromagProps import wavenumber
from Extractables.ElectromagProps import intrin_imp
from Extractables.ElectromagProps import appar_imp
from Extractables.ElectromagProps import appar_resis


miti = const.eps0.value
meab = const.mu0.value

angfreq = 0.00001
mu = 1
epsilon = 1
elecon = 1
prevappar = 1
thick=1

magperm = mu*meab

thick1 = 4.5
cond1 = 10**6
dieperm1 = 10*miti

wave1 = wavenumber(angfreq, magperm, cond1, dieperm1)
intrin1 = intrin_imp(angfreq, magperm, cond1, dieperm1)
apparimp1 = appar_imp(wave1, intrin1, intrin1, thick1)

print(wave1)
print(apparimp1)
print(intrin1)

print(np.tanh(5j))