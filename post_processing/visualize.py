
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import re

def create_animation(zeta_history, filename="vorticity_evolution.mp4", interval=50):
    fig, ax = plt.subplots()
    im = ax.imshow(zeta_history[0], origin='lower', cmap='seismic',
                   extent=[0, 2*np.pi, 0, 2*np.pi], animated=True)
    cbar = plt.colorbar(im, ax=ax)

    def update(frame):
        im.set_array(zeta_history[frame])
        ax.set_title(f"frame = {frame}")
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
    [f for f in os.listdir(output_dir) if f.endswith(".csv")],
    key=extract_step
)
zeta_history = []

for fname in files:
    data = np.loadtxt(os.path.join(output_dir, fname), delimiter=",")
    zeta_history.append(data)

create_animation(zeta_history, filename="vorticity_evolution.mp4", interval=50)
print("Video saved as vorticity_evolution.mp4")
