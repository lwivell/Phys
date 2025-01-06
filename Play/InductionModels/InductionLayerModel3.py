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

"""
Aim of this model, take arrays of depths for ocean and crusts, run those combinations, and plot the outcome as a time series, simulating a low flyby. 
So basically, initially, plot multiple on same plot, but get the maths working for taking the arrays of depths,
then worry about the plots over time. 

Need to switch from using lists to arrays, lists will only work for appending results, not input parameters. 

Build slow and test at each level so the dimensions of my maths are well understood.
"""