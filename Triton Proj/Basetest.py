import numpy as np 
import matplotlib.pyplot as plt
from Extractables.InductionModel import LayeredSystem
import astropy.constants as const

base = LayeredSystem(json=r"./Json Models/Triton.json")
big = LayeredSystem(json=r"./Json Models/Triton.json")
small = LayeredSystem(json=r"./Json Models/Triton.json")
miti = const.eps0.value
angfreqs = np.array(np.logspace(-9, 7, 800))
freqs = angfreqs/(2*np.pi)

low = 14.2
high = 14.9
lowp = low *60*60
highp = high*60*60
lowf = 1/lowp
highf = 1/highp
lowerf = 1/(141*60*60)

baseconds, j,k,l,m = base.iterate(angfreqs)
params = [135.9,150,160,175,192.5]
res = []
#for n in params:
#    base._structure[1][2] = n
#    base._structure[1][1] = n* 1000
#    Wl._structure[4][1] = (10.5 -(0.5*n)) * 1000
#    base._structure[2][1] = (334.4 - n) * 1000
#    a,b,c,d,e = base.iterate(angfreqs)
#    res.append(a)
big._structure[1][1] = 192.5 * 1000
big._structure[2][1] = 143.9 * 1000
big._structure[1][2] = 8

small._structure[1][1] = 135.9 * 1000
small._structure[2][1] = 200.5 * 1000
small._structure[1][2] = 1

a,s,d,f,g = big.iterate(angfreqs)
z,x,c,v,b = small.iterate(angfreqs)


fig, axs = plt.subplots(1,1)
axs.plot(freqs,baseconds, label = 'Base Parameters')
axs.plot(freqs,a,label='Biggest Signature')
axs.plot(freqs,z,label='Smallest Signature')
#axs.plot(freqs, res[0], label=params[0])
#axs.plot(freqs, res[1], label=params[1])
#axs.plot(freqs, res[2], label=params[2])
#axs.plot(freqs, res[3], label=params[3])
#axs.plot(freqs, res[4], label=params[4])
axs.fill_betweenx(
    y=np.linspace(min(baseconds), max(baseconds), 100),
    x1=lowf,
    x2=highf,
    color='blue',
    alpha=0.5,
    label="Shaded Range"
)
axs.set_xscale('log')
axs.set_yscale('log')
axs.vlines(x = lowerf, ymin=min(baseconds), ymax=max(baseconds*10))
# Set labels and limits
axs.set_xlabel('Frequency [Hz]')
axs.set_ylabel('Apparent Conductivity [S/m]')
#axs.set_xlim(10**(-3), 10)  # Uncomment and adjust as needed
#axs.set_ylim(10**(-4.5), 0.5)  # Uncomment and adjust as needed

# Add legend and grid
axs.legend()
axs.grid(True)

# Show the plot
plt.tight_layout()
plt.show()
