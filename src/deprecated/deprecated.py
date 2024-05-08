"""
Code for other approaches for solving the 1D TISE for bound states:
1. Numerov method forward integrator with 'wag the tail' energy search
2. scipy.integrate.odeint method with forward and backward integration to turning points
   followed by matching condition of logarithmic derivatives for energy search
"""

from scipy.integrate import odeint, simps, simpson
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# example potential functions
ho = lambda x: 0.5 * x ** 2
isw = lambda x: x * 0
lp = lambda x: np.heaviside(x, 1) * x


# numerov integrator
def integrator(lb, rb, energy, pot_func):
    x, dx = np.linspace(lb, rb, 500, retstep=True)
    phi = np.zeros_like(x)
    phi[0] = 0
    phi[1] = 0.0001
    g = 2 * (energy - pot_func(x))
    f = 1 + ((dx ** 2) / 12) * g
    for i in range(len(phi))[1:-1]:
        phi[i + 1] = ((12 - 10 * f[i]) * phi[i] - f[i - 1] * phi[i - 1]) / f[i + 1]
    # normalize with Simpsons
    norm_factor = simpson(abs(phi) ** 2, x)
    if abs(norm_factor) < 1e-8:
        print("Normalization factor is zero, setting to 1")
        norm_factor = 1
    return x, phi / np.sqrt(norm_factor)


def get_system_state(state, x, E, V):
    """
    Get the state vector at value of x in order to solve the
    system of first-order ODEs using scipy.integrate.odeint
    """

    y, z = state

    return np.array([z, 2 * (V(x) - E) * y])


def solve_se(x_vals, E, V):
    """
    Solve the 1D Time-independent Schrodinger Equation
        f''(x) = 2(V(x) - E)f(x)

    Where f(x) is the quantum wavefunction (AKA the 'state').
    To use scipy.integrate.odeint, turn this equation into a
    system of first order equations by defining

        y'(x) = z(x)
        z'(x) = g(x)y(x)

    where
        g(x) = 2(V(x) - E)

    For initial conditions (at xs[0]), we need y = 0 and z to
    be some small positive number, since the wavefunction decays
    to zero as x goes to plus/minus infinity
    """

    # find the classical turning points
    obj = V(x_vals) - E
    idx = [ind[0] for ind in np.argwhere(np.diff(np.sign(obj))) if abs(obj[ind]) > 1e-6]
    # print(f"Turning points are: {x_vals[idx]}")
    # calculate forward and backward integrated solutions
    psi_fwd = odeint(get_system_state, [0.0, 0.1], x_vals[:max(idx)], args=(E, V))[:, 0]
    psi_bwd = odeint(get_system_state, [0.0, -0.1], x_vals[min(idx):][::-1], args=(E, V))[:, 0][::-1]
    # count the number of nodes in forwardly integrated solution to fix sign in derivative
    node_idx = np.array([ind[0] for ind in np.argwhere(np.diff(np.sign(psi_fwd))) if abs(psi_fwd[ind]) > 1e-6])
    node_num = len(node_idx)
    psi_fwd *= (-1) ** node_num
    # concatenate forward and backward solutions together and return matching index
    sol = np.concatenate([psi_fwd[:min(idx)], psi_bwd])

    # normalization
    sol /= np.sqrt(simps(sol ** 2, x=x_vals))

    # compute the logarithmic derivative at matching point
    dx = x_vals[1] - x_vals[0]
    match_index = min(idx)
    gamma_left = (sol[match_index + 1] - sol[match_index - 1]) / (2 * dx * sol[match_index])
    gamma_right = (sol[match_index + 2] - sol[match_index]) / (2 * dx * sol[match_index + 1])
    delta = gamma_right - gamma_left

    return sol, delta, node_num


def find_energy_eigenvalue(E_l, n):
    delta_E = 0.1
    tolerance = 0.01
    x = np.linspace(-20, 20, int(5e3))
    dx = x[1] - x[0]
    pot_func = ho

    bracket_iterations = 0
    bisection_iterations = 0
    E_guess = E_l
    sol, diff, _ = solve_se(x, E=E_guess, V=pot_func)
    start_sign = np.sign(diff)
    if abs(diff) > tolerance:
        current_sign = start_sign
        # try the bracket the eigen-energy within the interval
        while current_sign == start_sign:
            E_l = E_guess
            E_guess += delta_E
            sol, diff, _ = solve_se(x, E=E_guess, V=pot_func)
            current_sign = np.sign(diff)
            bracket_iterations += 1
            print(f"Bracket iteration: {bracket_iterations}, Error {diff}")
            if current_sign != start_sign:
                E_u = E_guess
                break
        print(f"Found bracket. Lower: {E_l}, Upper: {E_u}")
        # start the bisection method
        while abs(diff) > tolerance:
            bisection_iterations += 1
            E_guess = (E_l + E_u) / 2
            sol, diff, nodes = solve_se(x, E_guess, V=pot_func)
            current_sign = np.sign(diff)
            if current_sign == -1:
                E_l = E_guess
            else:  # sign is 1
                E_u = E_guess
            print(f"Bisect iteration: {bisection_iterations}, Error {diff}")

    print(f"Found energy eigenvalue: {E_guess}, for {nodes}-th excited state")

    return E_guess, sol
