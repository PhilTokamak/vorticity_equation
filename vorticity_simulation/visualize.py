
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import glob

output_dir = "output"
files = sorted([f for f in os.listdir(output_dir) if f.endswith(".csv")])
frames = []

for fname in files:
    data = np.loadtxt(os.path.join(output_dir, fname), delimiter=",")
    fig, ax = plt.subplots(figsize=(6.08, 6.08))
    im = ax.imshow(data, cmap="RdBu", origin="lower", extent=[0, 2*np.pi, 0, 2*np.pi])
    ax.set_title(fname)
    plt.colorbar(im)
    plt.tight_layout()
    frame_name = f"frame_{fname.split('_')[1].split('.')[0]}.png"
    plt.savefig(frame_name)
    frames.append(imageio.v2.imread(frame_name))
    plt.close()


for filename in glob.glob("frame_*.png"):
    os.remove(filename)

imageio.mimsave("vorticity_evolution.mp4", frames, fps=10)
print("Video saved as vorticity_evolution.mp4")
