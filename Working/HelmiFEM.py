import skfem
from skfem.helpers import dot, grad  # helpers make forms look nice
@skfem.BilinearForm
def a(u, v, _):
    return dot(grad(u), grad(v))
import numpy as np
from helmi import Helmholtz
import matplotlib.pyplot as plt

# create a rectangular mesh with skfem:
x_pts = np.linspace(0, 200, 101)
y_pts = np.linspace(-100, 100, 101)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts)
mesh = mesh.with_subdomains({'air': lambda x: x[0] < 100,
                             'plastic': lambda x: x[0] >= 100})
mesh = mesh.with_boundaries({'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
                             'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
                             'bound_ymin': lambda x: np.isclose(x[1], y_pts[0]),
                             'bound_ymax': lambda x: np.isclose(x[1], y_pts[-1])})

element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

k0 = 0.3
eps_air = 1
mu_air = 1
eps_plastic = 2 - 0.1j
mu_plastic = 1
fem.assemble_subdomains(alpha={'air': 1 / mu_air, 
                               'plastic': 1 / mu_plastic}, 
                        beta={'air': -1 * k0 ** 2 * eps_air, 
                              'plastic': -1 * k0 ** 2 * eps_plastic}, 
                        f={'air': 0, 
                           'plastic': 0})

fem.assemble_boundaries_dirichlet(value={'bound_ymin':  0, 
                                         'bound_ymax':  0})

fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1 / mu_plastic * 1j * k0, 
                                   'bound_xmax': 1 / mu_plastic * 1j * k0}, 
                            q={'bound_xmin': 1 / mu_plastic * 2j * k0, 
                               'bound_xmax': 0})

fem.solve()
x_bound_xmin, y_bound_xmin = fem.basis.doflocs[:, fem.basis.get_dofs('bound_xmin')]

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as mplt

mag = np.sqrt((fem.phi_re**2) +(fem.phi_im**2))

fig, ax = mplt.subplots(1, 1)
#plot(fem.basis, fem.phi_re, ax=ax[0], colorbar=True)
#plot(fem.basis, fem.phi_im, ax=ax[1], colorbar=True)
plot(fem.basis, mag, colorbar=True)
#ax[0].set_aspect(1)
#ax[1].set_aspect(1)
mplt.tight_layout()
mplt.title('Air Plastic Example')
mplt.show()