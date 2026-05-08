import numpy as np
import matplotlib.pyplot as plt
import os

result_dir = "result"

Ns = [16, 32, 64, 128]
Nu_final = []

for N in Ns:
    filename = os.path.join(result_dir, f"nusselt_Ra_1e+04_N_{N}.txt")

    data = np.loadtxt(filename)

    # columns: step, time, Nu
    Nu = data[:, 2]

    # use last recorded value
    Nu_final.append(Nu[-1])

Nu_benchmark = 4.884

plt.figure(figsize=(7, 5))
plt.plot(Ns, Nu_final, "o-", linewidth=2, markersize=7, label="Computed")
plt.axhline(Nu_benchmark, linestyle="--", linewidth=2, label="Blankenbach Nu = 4.884")

plt.xlabel("Mesh size N")
plt.ylabel("Final Nusselt number")
plt.title(r"Mesh convergence for $Ra=10^4$")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(result_dir, "Nu_mesh_convergence_Ra_1e4.png"), dpi=300)
plt.show()

for N, Nu in zip(Ns, Nu_final):
    error = abs(Nu - Nu_benchmark)
    rel_error = error / Nu_benchmark * 100
    print(f"N={N:3d}, Nu={Nu:.6f}, error={error:.6f}, rel_error={rel_error:.2f}%")