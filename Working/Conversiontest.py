import numpy as np
from Extractables.Conversions import cartdepth_from_sph
from Extractables.Conversions import cartcon_from_sphcon

a = 1560000
r = (0.288,0.929, 0.994, 1)
cons = (10**6, 10**(-3), 1, (10**(-4)))

depths = []
cartcons =[]
for n in r:
    z = cartdepth_from_sph(a,n)
    depths.append(z)
    
    for m in cons:
        if r[x] == cons[m]:
            cartcon = cartcon_from_sphcon(m, n)
            cartcons.append(cartcon)

print(cartcons)
