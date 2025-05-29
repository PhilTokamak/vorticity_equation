import numpy as np
import matplotlib.pyplot as plt
import os
import re
from multigrid import *

def kinetic_energy_spectrum(ux, uy, dx):
    N = ux.shape[0]
    ux_hat = np.fft.fft2(ux)
    uy_hat = np.fft.fft2(uy)

    # Compute wave number grid
    kx = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    kx, ky = np.meshgrid(kx, ky, indexing='ij')
    k_mag = np.sqrt(kx**2 + ky**2)

    # Spectral energy density
    energy_density = 0.5 * (np.abs(ux_hat)**2 + np.abs(uy_hat)**2) / (N**2)

    # Flatten arrays
    k_flat = k_mag.ravel()
    e_flat = energy_density.ravel()

    # Shell average to get E(k)
    k_max = int(np.max(k_mag))
    k_bins = np.arange(0.5, k_max + 0.5, 1.0)  # Bin edges: [0.5, 1.5, 2.5, ...]
    k_vals = 0.5 * (k_bins[:-1] + k_bins[1:])

    # Histogram: sum of energies and counts
    E_sum, _ = np.histogram(k_flat, bins=k_bins, weights=e_flat)
    counts, _ = np.histogram(k_flat, bins=k_bins)

    # Average energy per shell
    E_k = np.zeros_like(counts, dtype=np.float64)
    nonzero = counts > 0
    E_k[nonzero] = E_sum[nonzero] / counts[nonzero]

    return k_vals, E_k


N = 128
L = 2 * np.pi
dx = L / N

output_dir = "output_Taylor-Green_fd"
def extract_step(filename):
    match = re.search(r"zeta_(\d+)\.csv", filename)
    return int(match.group(1)) if match else -1

files = sorted(
    [f for f in os.listdir(output_dir) if (f.endswith(".csv") and f.startswith("zeta"))],
    key=extract_step
)
zeta_history = []
psi_history = []
ux_history = []
uy_history = []

for fname in files:
    zeta = np.loadtxt(os.path.join(output_dir, fname), delimiter=",")
    psi = solve_poisson(zeta, dx)
    ux = - np.gradient(psi, axis=0) / dx
    uy = np.gradient(psi, axis=1) / dx
    zeta_history.append(zeta)
    psi_history.append(psi)
    ux_history.append(ux)
    uy_history.append(uy)

T = np.loadtxt(os.path.join(output_dir, "T.csv"), delimiter=",")

# Get k array
k_vals, _ = kinetic_energy_spectrum(ux_history[0], uy_history[0], dx)

E_k_history = []
for i in range(len(T)):
    _, E_k = kinetic_energy_spectrum(ux_history[i], uy_history[i], dx)
    E_k_history.append(E_k)
E_k_history = np.array(E_k_history)
k_index_to_plot = [0]
for ind in k_index_to_plot:
    plt.plot(T, E_k_history[:,ind], label = "E(k={})".format(k_vals[ind]))
plt.legend()
plt.savefig(output_dir + "/visualization/energy_evolution.pdf")
#plt.show()
plt.close()


#Time index to plot energy spectrum
t_ind_to_plot = [0, 230, 400]

for ind in t_ind_to_plot:
    plt.loglog(k_vals, E_k_history[ind,:], label='E(k) at t = {}'.format(T[ind]))
plt.loglog(k_vals, k_vals**(-5/3), '--', label=r'$k^{-5/3}$')  # Inverse cascade slope
plt.loglog(k_vals, k_vals**(-3), '--', label=r'$k^{-3}$')      # Forward enstrophy cascade
plt.legend()
plt.xlabel("Wavenumber k")
plt.ylabel("E(k)")
plt.title("Energy epectrum")
plt.savefig(output_dir + "/visualization/energy_spectrum.pdf")
#plt.show()






