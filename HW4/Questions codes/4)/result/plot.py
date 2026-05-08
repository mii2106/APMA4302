import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

result_dir = "result"

files = sorted(glob.glob(os.path.join(result_dir, "nusselt_Ra_*_N_*.txt")))

plt.figure(figsize=(8, 5))

for file in files:
    data = np.loadtxt(file)

    step = data[:, 0]
    time = data[:, 1]
    Nu = data[:, 2]

    name = os.path.basename(file)

    match = re.search(r"nusselt_Ra_(.*?)_N_(\d+)", name)
    if match:
        Ra = match.group(1)
        N = match.group(2)
        label = f"Ra={Ra}, N={N}"
    else:
        label = name

    plt.semilogx(time, Nu, linewidth=2, label=label)

plt.axhline(1.0, linestyle="--", linewidth=1, label="Nu = 1")
plt.xlabel("time")
plt.ylabel("Nusselt number")
plt.title("Nusselt number vs time")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(result_dir, "Nu_vs_time_all.png"), dpi=300)
plt.show()