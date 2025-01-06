import numpy as np
import matplotlib.pyplot as plt

def solve_wave_equation(c, dx, dt, L, T):
    """
    Solves the one-dimensional wave equation using the finite difference scheme.

    Parameters:
    - c: float
        Wave speed.
    - dx: float
        Spatial step size.
    - dt: float
        Time step size.
    - L: float
        Length of the domain.
    - T: float
        Total time.

    Returns:
    - numpy.ndarray:
        Array containing the solution of the wave equation at each time step.

    Raises:
    - ValueError:
        Raises an error if the stability condition is not satisfied.
    """

    # Calculate the number of spatial and temporal steps
    Nx = int(L / dx) + 1
    Nt = int(T / dt) + 1

    # Check the stability condition
    if c * dt / dx > 1:
        raise ValueError("Stability condition not satisfied. Reduce the time step size or increase the spatial step size.")

    # Initialize the solution array
    u = np.zeros((Nt, Nx))

    # Set the initial conditions
    x = np.linspace(0, L, Nx)
    u[0] = np.sin(2 * np.pi * x / L)

    # Solve the wave equation using the finite difference scheme
    for n in range(1, Nt):
        for i in range(1, Nx - 1):
            u[n, i] = 2 * (1 - (c * dt / dx) ** 2) * u[n - 1, i] - u[n - 2, i] + (c * dt / dx) ** 2 * (u[n - 1, i + 1] + u[n - 1, i - 1])

    return u

# Example usage of the solve_wave_equation function

# Parameters
c = 1.0  # Wave speed
dx = .05  # Spatial step size
dt = 0.00001  # Time step size
L = 1.0  # Length of the domain
T = 1.0  # Total time

# Solve the wave equation
solution = solve_wave_equation(c, dx, dt, L, T)

# Plot the solution
x = np.linspace(0, L, solution.shape[1])
t = np.linspace(0, T, solution.shape[0])
X, T = np.meshgrid(x, t)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, solution, cmap='viridis')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
plt.show()