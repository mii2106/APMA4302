import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

EXEC = "./reaction2d"

refine_levels = [2, 3, 4, 5, 6]
nprocs_list   = [1, 2, 4]

GAMMA = 100.0
P     = 3
ERR_RE   = re.compile(r"rel_error\s+\|u-uexact\|_2/\|uexact\|_2\s*=\s*([0-9.eE+\-]+)")
TIME_RE  = re.compile(r"^Time \(sec\):\s+([0-9.eE+\-]+)", re.MULTILINE)
SNES_RE  = re.compile(r"^\s*(\d+)\s+SNES Function norm", re.MULTILINE)

def run_case(nprocs, refine):
    cmd = [
        "mpirun", "-n", str(nprocs),
        EXEC,
        "-rct_gamma",  str(GAMMA),
        "-rct_p",      str(P),
        "-da_refine",  str(refine),
        "-snes_atol",  "1e-10",
        "-ksp_type",   "preonly",
        "-pc_type",    "lu",
        "-snes_monitor",
        "-log_view",
    ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    out = result.stdout + result.stderr

    # relative error
    m = ERR_RE.search(out)
    if not m:
        raise RuntimeError(f"No error found in output:\n{out}")
    relerr = float(m.group(1))

    # wall time
    m = TIME_RE.search(out)
    if not m:
        raise RuntimeError(f"No time found in output:\n{out}")
    walltime = float(m.group(1))

    # number of Newton iterations (last SNES iteration number)
    newton_iters = len(SNES_RE.findall(out)) - 1  # subtract iteration 0

    # grid size from refine level: base 9x9, each refine doubles-1
    nx = (9 - 1) * (2 ** refine) + 1

    return relerr, walltime, newton_iters, nx

# Run all cases
results = {}  

for nprocs, refine in product(nprocs_list, refine_levels):
    print(f"Running: nprocs={nprocs}, da_refine={refine} ...", end=" ", flush=True)
    try:
        relerr, walltime, newton_iters, nx = run_case(nprocs, refine)
        results[(nprocs, refine)] = (relerr, walltime, newton_iters, nx)
        print(f"nx={nx}, newton={newton_iters}, time={walltime:.3e}, relerr={relerr:.3e}")
    except Exception as e:
        print(f"FAILED: {e}")
        results[(nprocs, refine)] = (None, None, None, None)

print("\n" + "="*75)
print(f"{'nprocs':>8} {'refine':>8} {'grid':>10} {'Newton its':>12} {'Time (s)':>12} {'Rel. error':>12}")
print("="*75)
for nprocs in nprocs_list:
    for refine in refine_levels:
        relerr, walltime, newton_iters, nx = results[(nprocs, refine)]
        if relerr is not None:
            print(f"{nprocs:>8} {refine:>8} {nx:>4}x{nx:<4} {newton_iters:>12} {walltime:>12.3e} {relerr:>12.3e}")
        else:
            print(f"{nprocs:>8} {refine:>8} {'N/A':>10} {'N/A':>12} {'N/A':>12} {'N/A':>12}")
print("="*75)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# --- Time vs grid size ---
ax = axes[0]
for nprocs in nprocs_list:
    nxs   = [results[(nprocs, r)][3] for r in refine_levels if results[(nprocs, r)][0] is not None]
    times = [results[(nprocs, r)][1] for r in refine_levels if results[(nprocs, r)][1] is not None]
    ax.loglog(nxs, times, 'o-', label=f"{nprocs} proc(s)")
ax.set_xlabel("Grid size $n$", fontsize=12)
ax.set_ylabel("Wall time (s)", fontsize=12)
ax.set_title("Wall time vs grid size", fontsize=12)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend()

# --- Relative error vs grid size ---
ax = axes[1]
for nprocs in nprocs_list[:1]:  # error doesn't depend on nprocs
    nxs    = [results[(nprocs, r)][3] for r in refine_levels if results[(nprocs, r)][0] is not None]
    errors = [results[(nprocs, r)][0] for r in refine_levels if results[(nprocs, r)][0] is not None]
    ax.loglog(nxs, errors, 'o-', color='steelblue', label="Rel. error")
    # fit slope
    log_n = np.log(nxs)
    log_e = np.log(errors)
    slope = np.polyfit(log_n, log_e, 1)[0]
    ax.loglog(nxs, errors[0] * (np.array(nxs)/nxs[0])**slope,
              '--', color='gray', label=f"slope = {slope:.2f}")
ax.set_xlabel("Grid size $n$", fontsize=12)
ax.set_ylabel("Relative error", fontsize=12)
ax.set_title("Convergence vs grid size", fontsize=12)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend()

# --- Newton iterations vs grid size ---
ax = axes[2]
for nprocs in nprocs_list[:1]:
    nxs     = [results[(nprocs, r)][3]    for r in refine_levels if results[(nprocs, r)][0] is not None]
    newtons = [results[(nprocs, r)][2]    for r in refine_levels if results[(nprocs, r)][2] is not None]
    ax.plot(nxs, newtons, 'o-', color='darkorange', label="Newton its")
ax.set_xlabel("Grid size $n$", fontsize=12)
ax.set_ylabel("Newton iterations", fontsize=12)
ax.set_title("Newton iterations vs grid size", fontsize=12)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend()

plt.suptitle(r"Scaling study: $\gamma=100$, $p=3$", fontsize=13)
plt.tight_layout()
plt.savefig("scaling_study.png", dpi=150)
plt.show()
print("Plot saved as scaling_study.png")