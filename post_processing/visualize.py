
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import re

def create_animation(zeta_history, T, filename="vorticity_evolution.mp4", interval=50):
    fig, ax = plt.subplots()
    im = ax.imshow(zeta_history[0], origin='lower', cmap='seismic',
                   extent=[0, 2*np.pi, 0, 2*np.pi], animated=True)
    cbar = plt.colorbar(im, ax=ax)

    def update(frame):
        im.set_array(zeta_history[frame])
        ax.set_title(f"t = {T[frame]}")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=len(zeta_history),
                                  interval=interval, blit=True)

    ani.save(filename, writer='ffmpeg', dpi=150)
    plt.close()


output_dir = "output"
def extract_step(filename):
    match = re.search(r"zeta_(\d+)\.csv", filename)
    return int(match.group(1)) if match else -1

files = sorted(
    [f for f in os.listdir(output_dir) if (f.endswith(".csv") and f.startswith("zeta"))],
    key=extract_step
)
zeta_history = []

for fname in files:
    data = np.loadtxt(os.path.join(output_dir, fname), delimiter=",")
    zeta_history.append(data)

T = np.loadtxt(os.path.join(output_dir, "T.csv"), delimiter=",")

create_animation(zeta_history, T, filename="vorticity_evolution.mp4", interval=50)
print("Video saved as vorticity_evolution.mp4")


# plot mean vorticity evolution
diag_filename = "mean_vorticity.csv"
zeta_mean = np.loadtxt(os.path.join(output_dir, "diagnostics", diag_filename), delimiter=",")

fig, ax = plt.subplots()
ax.plot(T, zeta_mean)
ax.set_ylabel(r"Mean $\zeta$")
ax.set_xlabel("Time")
ax.set_title(r"Mean Vorticity Evolution")
ax.grid()
plt.savefig("mean_vorticity.pdf")
plt.close()

# plot mean kinetic energy evolution
diag_filename = "mean_kinetic_E.csv"
K_mean = np.loadtxt(os.path.join(output_dir, "diagnostics", diag_filename), delimiter=",")

fig, ax = plt.subplots()
ax.plot(T, K_mean)
ax.set_ylabel(r"Mean E_k")
ax.set_xlabel("Time")
ax.set_title(r"Mean Kinetic Energy Evolution")
ax.grid()
plt.savefig("mean_kinetic_E.pdf")
plt.close()

# plot mean enstrophy evolution
diag_filename = "mean_enstrophy.csv"
enst_mean = np.loadtxt(os.path.join(output_dir, "diagnostics", diag_filename), delimiter=",")

fig, ax = plt.subplots()
ax.plot(T, enst_mean)
ax.set_ylabel(r"Mean Enstrophy")
ax.set_xlabel("Time")
ax.set_title(r"Mean Enstrophy Evolution")
ax.grid()
plt.savefig("mean_enstrophy.pdf")
plt.close()
