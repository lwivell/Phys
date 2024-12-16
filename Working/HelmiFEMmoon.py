import skfem
from skfem.helpers import dot, grad  # helpers make forms look nice
@skfem.BilinearForm
def a(u, v, _):
    return dot(grad(u), grad(v))
import numpy as np
from helmi import Helmholtz
import matplotlib.pyplot as plt
import astropy.constants as const

miti = const.eps0.value
meab = const.mu0.value

# create a rectangular mesh with skfem:
x_pts = np.linspace(0, 100, 151)
y_pts = np.linspace(0, 10, 51)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts)
mesh = mesh.with_subdomains({'ice': lambda x: x[0] < 26,
                             'water': lambda x: x[0] >= 26})
mesh = mesh.with_boundaries({'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
                             'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
                             'bound_ymin': lambda x: np.isclose(x[1], y_pts[0]),
                             'bound_ymax': lambda x: np.isclose(x[1], y_pts[-1])})

element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

angfreq = 5000

eps_ice = 3.5*miti
mu_ice = 1*meab
cond_ice = 0.0001
k_ice = 1j*mu_ice*angfreq * (cond_ice + (1j*angfreq*eps_ice))/0.001

eps_water = 85*miti
mu_water = 1*meab
cond_water = 1
k_water = 1j*mu_water*angfreq * (cond_water + (1j*angfreq*eps_ice))/0.001

eps_rock = 5*miti
mu_rock = 1*meab
cond_rock = 0.001
k_rock = 1j*mu_rock*angfreq * (cond_rock + (1j*angfreq*eps_rock))/0.001

fem.assemble_subdomains(alpha={'ice': 1, 
                               'water': 1}, 
                        beta={'ice': -1 * k_ice, 
                              'water': -1 * k_water}, 
                        f={'ice': 0, 
                           'water': 0})

#fem.assemble_boundaries_dirichlet(value={'bound_ymin':  1, 
#                                         'bound_ymax':  1})

fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1j *np.sqrt(k_ice), 
                                   'bound_xmax': 1j * np.sqrt(k_water)}, 
                            q={'bound_xmin': 2j * (angfreq/(3e8 * 0.001)), 
                               'bound_xmax': -1j * (angfreq/(3e8 * 0.001))})

fem.solve()
x_bound_xmin, y_bound_xmin = fem.basis.doflocs[:, fem.basis.get_dofs('bound_xmin')]

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as mplt

mag = np.sqrt((fem.phi_re**2) +(fem.phi_im**2))

fig, ax = mplt.subplots(3, 1)
plot(fem.basis, fem.phi_re, ax=ax[0], colorbar=True)
plot(fem.basis, fem.phi_im, ax=ax[1], colorbar=True)
plot(fem.basis, mag, ax=ax[2], colorbar=True)
ax[0].set_xlim(0,61)
ax[1].set_xlim(0,61)
ax[2].set_xlim(0,61)
#ax[0].set_aspect(1)
#ax[1].set_aspect(1)
mplt.tight_layout()


#mplt.title('Air Plastic Example')
mplt.show()