import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

miti = const.eps0.value
meab = const.mu0.value

omega = 0
mu = 1
epsilon = 1
cond = 1
prevappar = 1
h=1

k = (((omega**2)*mu*meab*epsilon*miti)-((1j)*omega*mu*meab*cond))**(1/2)
lowk = (-(1j)*omega*mu*meab*cond)**(1/2)
highk = ((omega**2)*mu*meab*epsilon*miti)


intrinimp = (((1j)*mu*meab*omega)/(cond+((1j)*omega*epsilon*miti)))**(1/2)

apparimp = intrinimp*((prevappar+(intrinimp*np.tanh((1j)*k*h)))/(intrinimp+(prevappar*np.tanh((1j)*k*h))))

