import matplotlib.pyplot as plt

iters = list(range(9))
res = [
    2.097782055445e+03,
    1.855374980628e+03,
    1.627399462002e+03,
    1.254451498513e+03,
    9.080764440243e+02,
    9.415433107178e+01,
    1.442912640207e+00,
    3.570373316336e-04,
    3.917317745967e-11
]

plt.semilogy(iters, res, marker='o')
plt.xlabel('Newton iteration')
plt.ylabel(r'$\|F(u)\|_2$')
plt.title(r'Residual $\gamma=100$, $p=3$')
plt.grid(True)
plt.show()