import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt

EXEC = "./bvp"  

ms = [40, 80, 160, 320, 640, 1280]
ks = [1, 5, 10]
gamma = 0.0
c = 3.0

ERR_RE = re.compile(r"system is =\s*([0-9.+-eE]+)")
def run_case(m, k):
    cmd = [
        EXEC,
        f"-bvp_m", str(m),
        f"-bvp_gamma", str(gamma),
        f"-bvp_k", str(k),
        f"-bvp_c", str(c),
    ]
    out = subprocess.check_output(cmd, text=True)
    match = ERR_RE.search(out)
    if not match:
        raise RuntimeError(f"No pude encontrar error relativo en la salida:\n{out}")
    relerr = float(match.group(1))
    h = 1.0/(m-1)  
    return h, relerr

def fit_order(hs, errs):
    x = np.log(hs)
    y = np.log(errs)
    p = np.polyfit(x, y, 1)  
    order = p[0]
    return order, p

results = {}

for k in ks:
    hs = []
    errs = []
    for m in ms:
        h, e = run_case(m, k)
        hs.append(h)
        errs.append(e)
        print(f"k={k:2d}, m={m:4d}, h={h:.3e}, relerr={e:.3e}")
    hs = np.array(hs)
    errs = np.array(errs)
    order, p = fit_order(hs, errs)
    results[k] = (hs, errs, order)

# plot log-log
plt.figure(figsize=(8,6))
for k in ks:
    hs, errs, order = results[k]
    plt.loglog(hs, errs, marker='o', label=f"k={k}, slope = {order:.2f}")

plt.xlabel("h")
plt.ylabel("relative error")
plt.title("Convergence of error vs h (gamma=0)")
plt.grid(True, which="both")
plt.legend()
plt.show()