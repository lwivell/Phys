import numpy as np
import matplotlib.pyplot as plt

def magnetic_dipole_field(x, y, m=7.65e15):
    """
    Calculate the magnetic field components (Bx, By) for a magnetic dipole at the origin.
    :param x: x-coordinate(s) where the field is calculated
    :param y: y-coordinate(s) where the field is calculated
    :param m: Dipole moment (default is 1)
    :return: Tuple of arrays (Bx, By)
    """
    r = np.sqrt(x**2 + y**2)
    r_cubed = r**3
    with np.errstate(divide='ignore', invalid='ignore'):
        Bx = (3 * x * y * m) / r_cubed
        By = (m * (3 * y**2 - r**2)) / r_cubed
        Bx[np.isnan(Bx)] = 0  # Handle NaN at the origin
        By[np.isnan(By)] = 0
    return Bx, By
r_moon = 1.737e6

# Create a grid of points
x = np.linspace(-2*r_moon, 2*r_moon, 100)
y = np.linspace(-2*r_moon, 2*r_moon, 100)
x, y = np.meshgrid(x, y)

# Calculate the field components
Bx, By = magnetic_dipole_field(x, y)

# Plot the field lines
plt.figure(figsize=(8, 8))
plt.streamplot(x, y, Bx, By, color=np.sqrt(Bx**2 + By**2), linewidth=1.5, cmap='plasma')
plt.scatter(0, 0, color='red', s=100, label='Magnetic Dipole')
plt.title('Magnetic Dipole Field')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.legend()
plt.colorbar(label='Field Strength')
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()
