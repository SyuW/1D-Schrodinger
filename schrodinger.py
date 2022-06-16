from scipy.integrate import odeint, simps, simpson
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# example potential functions
ho = lambda x: 0.5 * x ** 2
isw = lambda x: x * 0
lp = lambda x: np.heaviside(x, 1) * x


def find_eigenstates_energies(n_eigenstates_to_find, V):
    """
    Discretize the x-axis, then using a finite difference definition
    of derivative: f'(x) = (f(x+h) - f(x-h)) / 2h, we have a matrix
    representation D of the derivative operator acting on discretized
    eigenfunction. This gives us the eigenvalue problem to solve using
    available diagonalization routines:

    H * v = E * v

    where H is the Hamiltonian defined by

    H = -0.5 * D^2 + V(x)
    """

    num_points = 1000
    x_vals = np.linspace(-10, 10, num_points)

    dx = x_vals[1] - x_vals[0]
    off = np.ones(num_points - 1)  # off diagonals
    on = np.ones(num_points)  # on diagonal
    d_squared = 0.5 * (np.diag(off, k=1) - 2 * np.diag(on, k=0) + np.diag(off, k=1)) / (dx ** 2)

    # Hamiltonian
    H = -d_squared + np.diag(V(x_vals))

    # solve for eigenvalues/eigenstates and sort
    energies, eigenstates = np.linalg.eig(H)
    idx = energies.argsort()
    energies = energies[idx]
    eigenstates = eigenstates[:, idx]

    # normalize with simpson's rule
    for ii in range(n_eigenstates_to_find):
        eigenstates[:, ii] /= simps(eigenstates[:, ii] ** 2, x=x_vals)

    return energies[:n_eigenstates_to_find], eigenstates[:, n_eigenstates_to_find]


if __name__ == "__main__":
    # ----- example potential functions ----- #
    # harmonic oscillator
    ho = lambda x: 0.5 * x ** 2
    # square well
    isw = lambda x, w, U: U * (abs(x) > w)
    # linear potential well
    tp = lambda x, k: k * abs(x)
    # quartic oscillator
    qo = lambda x: 0.5 * x ** 4
    # anharmonic oscillator
    aho = lambda x: 0.25 * x ** 2 + 0.5 * x ** 4
    # linear ramp
    lr = lambda x, k: k * x * (x > 0)
    # ------------------------------- #

    parser = argparse.ArgumentParser(description="Enter a one-dimensional potential")
    parser.add_argument("--potential", metavar="0.5*x**2", help="potential")
    parser.add_argument("--N", metavar=10, help="Number of energy states to find")
    parser.add_argument("--range", help="Size of domain")
    parser.add_argument("-v", type=bool, metavar="N", help="verbose")

    potential_func = lambda x: eval(args.v)

    args = parser.parse_args()
