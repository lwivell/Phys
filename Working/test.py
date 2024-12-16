import skfem as fem
from skfem.helpers import dot, grad  # helpers make forms look nice
@fem.BilinearForm
def a(u, v, _):
    return dot(grad(u), grad(v))