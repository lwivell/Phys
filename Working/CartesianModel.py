import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from Extractables.ElectromagProps import wavenumber
from Extractables.ElectromagProps import intrin_imp
from Extractables.ElectromagProps import appar_imp
from Extractables.ElectromagProps import appar_resis


miti = const.eps0.value
meab = const.mu0.value

angfreq = 0
mu = 1
epsilon = 1
elecon = 1
prevappar = 1
thick=1

dieperm = epsilon*miti
magperm = mu*meab

