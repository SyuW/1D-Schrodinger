from scipy.integrate import solve_bvp
from scipy.integrate import simpson
from scipy.optimize import curve_fit, fsolve
from scipy.special import hermite
from math import factorial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Harmonic oscillator potential: p^2/2m + x^2/2
ho = lambda x: 0.5 * x ** 2

# numerov integrator
def integrator(lb, rb, energy, pot_func):
    x, dx = np.linspace(lb, rb, int(5e3), retstep=True)
    phi = np.zeros_like(x)
    phi[0] = 0
    phi[1] = 1e-4
    g = 2 * (energy - pot_func(x))
    f = 1 + g*dx**2/12
    
    for i in range(len(phi))[1:-1]:
        phi[i+1] = ((12 - 10 * f[i]) * phi[i] - f[i-1] * phi[i-1]) / f[i+1]
    
    # normalize with Simpsons
    norm_factor = simpson(abs(phi)**2, x)
    if abs(norm_factor) < 1e-8:
        print("Normalization factor is zero, setting to 1")
        norm_factor = 1
    return x, phi / np.sqrt(norm_factor)

# plot for nth excited state
n = 3
E = n + 0.5
# exact solution using Hermite polynomials
H_n = hermite(n=n)
x, phi = integrator(lb=-6, rb=6, energy=E, pot_func=ho)
exact_sol = ((np.pi)**(-0.25)) * (1 / np.sqrt(2 ** n * factorial(n)) * H_n(x) * np.exp(-x ** 2 / 2))

plt.title("Exact versus numerical comparison")
plt.plot(x, (-1) ** n * phi, color="cyan", label="numerical")
plt.plot(x, exact_sol, color="red", label="exact")
plt.legend()
plt.show()