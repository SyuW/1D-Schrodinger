from scipy.integrate import odeint, simps, simpson
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib.cm import get_cmap


def find_eigenstates_energies(n_eigenstates_to_find, V, range_size):
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
    x_vals = np.linspace(-range_size, range_size, num_points)

    dx = x_vals[1] - x_vals[0]
    off = np.ones(num_points - 1)  # off diagonals
    on = np.ones(num_points)  # on diagonal
    d_squared = 0.5 * (np.diag(off, k=1) - 2 * np.diag(on, k=0) + np.diag(off, k=-1)) / (dx ** 2)

    # Hamiltonian
    H = -d_squared + np.diag(V(x_vals))

    # solve for eigenvalues/eigenstates and sort
    energies, eigenstates = np.linalg.eig(H)
    idx = energies.argsort()
    energies = energies[idx]
    eigenstates = eigenstates[:, idx]

    # normalize with simpson's rule
    for ii in range(n_eigenstates_to_find):
        eigenstates[:, ii] /= np.sqrt(simps(eigenstates[:, ii] ** 2, x=x_vals))
        print(f"Energy Level {ii}: {energies[ii]}")

    return energies[:n_eigenstates_to_find], eigenstates[:, :n_eigenstates_to_find]


def draw_double_figure(energies, eigenstates, V, range_size):
    """
    Draw a double figure with left plot containing energy spectrum and right plot
    containing stationary states for given potential
    """

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=False)
    num_points = 1000
    xgrid = np.linspace(-range_size, range_size, num_points)

    # set the x-axis limits based on turning points
    top = np.max(energies) * 1.5
    ax1.set_title("Energy Spectrum")
    ax1.set_xlabel("Position, x")
    ax1.set_ylabel("Energy, E")
    ax1.set_ylim([0, top])
    ax1.set_xlim([-range_size, range_size])
    ax1.plot(xgrid, V(xgrid), color="red", zorder=5)   # plot the potential

    ax2.set_title("Probability densities of eigenstates")
    ax2.set_xlabel("Position, x")
    ax2.set_ylabel("Probability density")

    # find the classical turning points of the potential for topmost line
    spacing = np.max(energies) / len(energies)
    topmost = np.min(energies) + len(energies) * spacing
    idx = np.argwhere(np.diff(np.sign(V(xgrid)-topmost)))
    x_l = np.min(xgrid[idx])
    x_r = np.max(xgrid[idx])
    lb = x_l - abs(x_l - x_r) / 2
    rb = x_r + abs(x_l - x_r) / 2
    # set the x-axis limits based on turning points
    ax2.set_xlim([lb, rb])
    ax2.set_ylim([0, topmost+spacing])

    color_map = get_cmap("inferno")
    ax2.plot(xgrid, V(xgrid), color="red", zorder=5)   # plot the potential
    for ii, energy in enumerate(energies):
        ax1.axhline(y=energy, linestyle="--", color=color_map(energy / top), label=f"E{ii}")
        base = np.min(energies) + ii * spacing
        ax2.axhline(y=base, color="black", linestyle="--", zorder=15)
        ax2.plot(xgrid, base + eigenstates[:, ii] ** 2, color="blue", zorder=10)
        ax2.fill_between(xgrid, base, base + eigenstates[:, ii] ** 2, color="blue", alpha=0.4)
        ax2.annotate(rf"$|\psi_{{{ii}}}|^2$", (lb+0.5, base-spacing*0.25), size=7.5)

    ax1.legend(loc="upper left", prop={'size': 6})

    return fig


def draw_single_figure():
    """

    """

    return


def draw_superposition_animation(energies, eigenstates, coeffs):
    """

    """

    fig, ax = plt.subplots()

    return fig


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
    # morse potential
    mp = lambda x: 30 * (1 - np.exp(-1 * (x - 1))) ** 2
    # ------------------------------- #

    parser = argparse.ArgumentParser(description="Enter a one-dimensional potential")
    parser.add_argument("--potential", metavar="0.5*x**2", type=str, required=True,
                        help="Potential")
    parser.add_argument("-c", "--coeffs", nargs="+", type=float, required=True,
                        help="Superposition coefficients")
    parser.add_argument("-r", "--range", type=float, nargs="?", required=False, const=1, default=10,
                        help="Size of domain")
    parser.add_argument("-N", type=int, metavar=8, nargs="?", required=False, const=1, default=8,
                        help="Number of energy states to find")
    parser.add_argument("--double_fig", type=bool, nargs="?", required=False, const=1, default=True,
                        help="Create a double figure")
    parser.add_argument("-v", type=bool, metavar="N", nargs="?", required=False, const=1, default=True,
                        help="Verbose mode")
    # get the command line arguments
    args = parser.parse_args()

    coeffs = np.array(args.coeffs)
    coeffs /= np.linalg.norm(coeffs)
    pot_func = lambda x: eval(args.potential)

    print(coeffs)
    pass

    v, w = find_eigenstates_energies(n_eigenstates_to_find=4, V=pot_func, range_size=args.range)
    finished_plot = draw_double_figure(v, w, V=pot_func, range_size=args.range)

    finished_plot.show()

    pass