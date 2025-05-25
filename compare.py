import os
import numpy as np
import re

output_dir_py = "output_py"
output_dir_cpp = "vorticity_simulation/build/output"
def extract_step(filename):
    match = re.search(r"zeta_(\d+)\.csv", filename)
    return int(match.group(1)) if match else -1

files = sorted(
    [f for f in os.listdir(output_dir_py) if f.endswith(".csv")],
    key=extract_step
)
zeta_history = []

for fname in files:
    data_py = np.loadtxt(os.path.join(output_dir_py, fname), delimiter=",")
    data_cpp = np.loadtxt(os.path.join(output_dir_cpp, fname), delimiter=",")
    if np.max(abs(data_py - data_cpp)) < 1e-10:
        print(fname, "ok")
    else:
        print(fname, "different!",
              np.argmax(abs(data_py - data_cpp))//data_py.shape[0],
              np.argmax(abs(data_py - data_cpp))%data_py.shape[0],
              np.max(abs(data_py - data_cpp)))