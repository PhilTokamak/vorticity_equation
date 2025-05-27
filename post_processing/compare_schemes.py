import numpy as np
import matplotlib.pyplot as plt
import os

output_dirs = ["output_arakawa-2", "output_arakawa-4", "output_fd"]
scheme_name = ["Arakawa-2", "Arakawa-4", "Centered Difference"]

zeta_mean = []
K_mean = []
enst_mean = []

# Read time points
T = np.loadtxt(os.path.join(output_dirs[0], "T.csv"), delimiter=",")


# plot mean vorticity evolution
for output_dir in output_dirs:
    zeta_mean.append(np.loadtxt(
        os.path.join(output_dir, "diagnostics", "mean_vorticity.csv"), delimiter=","))
    K_mean.append(np.loadtxt(
        os.path.join(output_dir, "diagnostics", "mean_kinetic_E.csv"), delimiter=","))
    enst_mean.append(np.loadtxt(
        os.path.join(output_dir, "diagnostics", "mean_enstrophy.csv"), delimiter=","))

fig, ax = plt.subplots()
for i in range(len(zeta_mean)):
    ax.plot(T, zeta_mean[i], label=scheme_name[i])
ax.set_ylabel(r"$\zeta$")
ax.set_xlabel("time")
ax.set_title(r"Mean Vorticity Evolution")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("mean_vorticity_compare.pdf")
plt.close()

# plot mean kinetic energy evolution
fig, ax = plt.subplots()
for i in range(len(K_mean)):
    ax.plot(T, K_mean[i], label=scheme_name[i])
ax.set_ylabel(r"E_k")
ax.set_xlabel("time")
ax.set_title(r"Mean Kinetic Energy Evolution")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("mean_kinetic_E_compare.pdf")
plt.close()

# plot mean enstrophy evolution
fig, ax = plt.subplots()
for i in range(len(enst_mean)):
    ax.plot(T, enst_mean[i], label=scheme_name[i])
ax.set_ylabel(r"Enstrophy")
ax.set_xlabel("time")
ax.set_title(r"Mean Enstrophy Evolution")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("mean_enstrophy_compare.pdf")
plt.close()