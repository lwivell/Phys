import skfem
from skfem.helpers import dot, grad  # helpers make forms look nice
@skfem.BilinearForm
def a(u, v, _):
    return dot(grad(u), grad(v))
import numpy as np
from helmi import Helmholtz
import matplotlib.pyplot as plt

# create a rectangular mesh with skfem:
x_pts = np.linspace(0, 122, 123)
y_pts = np.linspace(-200, 200, 201)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts)
mesh = mesh.with_subdomains({'air': lambda x: x[0] < 52,
                             'plastic': lambda x: x[0] >= 52})
mesh = mesh.with_boundaries({'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
                             'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
                             'bound_ymin': lambda x: np.isclose(x[1], y_pts[0]),
                             'bound_ymax': lambda x: np.isclose(x[1], y_pts[-1])})

element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

f = 1000
meshsize = 500
k0 = ((2*np.pi*f)/(3e8 / meshsize))
eps_air = 3.5
mu_air = 1
eps_plastic = 85
mu_plastic = 1
fem.assemble_subdomains(alpha={'air': 1 / mu_air, 
                               'plastic': 1 / mu_plastic}, 
                        beta={'air': -1 * k0 ** 2 * eps_air, 
                              'plastic': -1 * k0 ** 2 * eps_plastic}, 
                        f={'air': 0, 
                           'plastic': 0})

fem.assemble_boundaries_dirichlet(value={'bound_ymin':  -0.1, 
                                         'bound_ymax':  -0.1})

fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1 / mu_plastic * 1j * k0, 
                                   'bound_xmax': 1 / mu_plastic * 1j * k0}, 
                            q={'bound_xmin': 1 / mu_plastic * 2j * k0, 
                               'bound_xmax': 0})

fem.solve()
x_bound_xmin, y_bound_xmin = fem.basis.doflocs[:, fem.basis.get_dofs('bound_xmin')]

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as mplt

mag = np.sqrt((fem.phi_re**2) +(fem.phi_im**2))

fig, ax = mplt.subplots(1, 2)
plot(fem.basis, fem.phi_re,ax=ax[0], colorbar=True)
plot(fem.basis, fem.phi_im, ax=ax[1], colorbar=True)
#plot(fem.basis, mag, colorbar=True)
#ax[0].set_aspect(1)
#ax[1].set_aspect(1)
mplt.tight_layout()
mplt.title('Air Plastic Example')
mplt.show()