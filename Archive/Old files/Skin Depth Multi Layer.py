import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

# Constants
mu_0 = const.mu0.value # Permeability of free space in H/m

# Define layer properties: conductivity (S/m), thickness (m)
layers = [
    {'conductivity': 1e-4, 'thickness': 26000},  # Layer 1 (Top Layer)
    {'conductivity': 1, 'thickness': 35000},  # Layer 2
    {'conductivity': 1e-3, 'thickness': 190000}   # Layer 3 (Bottom Layer)
]

# Function to calculate skin depth for a single layer
def skin_depth(conductivity, frequency, mu=mu_0):
    omega = 2 * np.pi * frequency  # Angular frequency
    return np.sqrt(2 / (mu * conductivity * omega))

# Function to calculate effective skin depth for a three-layer medium
def effective_skin_depth(frequency, layers):
    skin_depths = [skin_depth(layer['conductivity'], frequency) for layer in layers]
    thicknesses = [layer['thickness'] for layer in layers]
    
    # Harmonic mean approach for effective skin depth
    effective_depth_inv = sum(d / (delta ** 2) for d, delta in zip(thicknesses, skin_depths))
    effective_skin_depth = 1 / np.sqrt(effective_depth_inv)
    return effective_skin_depth

# Frequency range for the plot (10 Hz to 10 MHz)
frequencies = np.logspace(-7, 4, 500)  # 10 Hz to 10 MHz
effective_depths = [effective_skin_depth(f, layers) for f in frequencies]

# Plotting
plt.figure(figsize=(10, 6))
plt.loglog(frequencies, effective_depths, label="Effective Skin Depth (3-layer medium)")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Effective Skin Depth [km]")
plt.title("Effective Skin Depth at Europa (Average Model)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.show()

datas = [frequencies, effective_depths]

np.save(r".\Enceladus Models\EffectiveSkinDepth.npy", datas)