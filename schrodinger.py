from scipy.integrate import odeint, simps, simpson
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


class Hamiltonian():

    def __init__(self, potential, domain):
        """
        Constructor for the Hamiltonian class

        potential should be a callable function representing the potential energy

        """
        self.V = potential
        self.domain = domain


    # solve method for solving the Hamiltonian
    def tise_solve(self, n_eigenstates_to_find):
        """
        Solve the time-independent version of the Schrodinger equation.
        Discretize the x-axis, then using a finite difference definition of the
        derivative: f'(x) = (f(x+h) - f(x-h)) / 2h, we have a matrix representation D
        of the derivative operator acting on a discretized eigenfunction, allowing us
        to represent the kinetic term in the Hamiltonian in the position basis. This
        gives us an eigenvalue problem to solve using available diagonalization routines:

        H * v = E * v

        where H is the Hamiltonian defined by

        H = -0.5 * D^2 + V(x)
        """

        num_points = 1000
        x_vals = np.linspace(self.domain[0], self.domain[1], num_points)

        dx = x_vals[1] - x_vals[0]
        off = np.ones(num_points - 1) # off diagonals
        on = np.ones(num_points) # on diagonal
        d_squared = 0.5 * (np.diag(off, k=1) - 2 * np.diag(on, k=0) + np.diag(off, k=-1)) / (dx ** 2)

        # Hamiltonian
        H = -d_squared + np.diag(self.V(x_vals))

        # solve for eigenvalues/eigenstates and sort in order of increasing energy
        energies, eigenstates = np.linalg.eig(H)
        idx = energies.argsort()
        energies = energies[idx]
        eigenstates = eigenstates[:, idx]

        # normalize with simpson's rule
        for ii in range(n_eigenstates_to_find):
            eigenstates[:, ii] /= np.sqrt(simps(eigenstates[:, ii] ** 2, x=x_vals))
            print(f"Energy Level {ii}: {energies[ii]}")

        return energies[:n_eigenstates_to_find], eigenstates[:, :n_eigenstates_to_find]


# helper function for drawing a state onto an Axes object when updating through key press
def plot_state(linedata, state):
    # note that there is a modulus squared because of the Born rule
    linedata


def draw_figure(energies, eigenstates, V, domain):
    """
    Draw a double figure with left plot containing the energy spectrum and right
    plot containing the stationary states for a given potential
    """

    num_eigenstates = len(energies)

    # state index for updating the current plotted eigenstate 
    # through a key press. Wrap in a dictionary to make it mutable (hackily).
    whats_the_current_state = {"current": 0}

    # for updating the plot
    def on_press(event):
        if event.key in ["left", "right"]:
            state_index = whats_the_current_state["current"]
            if event.key == 'left':
                inc = -1
            if event.key == 'right':
                inc = 1
            # first change the spectrum plot, meaning change the selected line's color to red 
            spectral_lines[state_index].set_color("black")
            state_index = (state_index + inc) % num_eigenstates
            spectral_lines[state_index].set_color("red")
            # second, change the eigenstates plot
            current_state.set_ydata(eigenstates[:, state_index] ** 2 + energies[state_index])
            fill = fill_dict["fill"]
            fill.remove()
            fill = ax2.fill_between(xgrid, energies[state_index],
                                    energies[state_index] + eigenstates[:, state_index] ** 2,
                                    color="blue", alpha=0.4)
            fill_dict["fill"] = fill
            whats_the_current_state["current"] = state_index
            ax2.set_title(rf"Probability density of eigenstate $|\psi_{{{state_index}}}|^2$")
            fig.canvas.draw()

    # set up the figure and axes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,8), 
                                   gridspec_kw={'width_ratios': [1,4]}, sharey=True)
    fig.canvas.mpl_connect('key_press_event', on_press)

    num_points = 1000
    xgrid = np.linspace(domain[0], domain[1], num_points)

    # spectrum plot
    top = np.max(energies) * 1.5
    ax1.set_ylabel("Energy, E")
    ax1.set_ylim([0, top])
    # don't need the x-axis ticks and labels
    ax1.set_xticks([])
    ax1.set_xticklabels([])

    spectral_lines = []
    for E in energies:
        spec = ax1.axhline(y=E, linestyle="--", color="black")
        spectral_lines.append(spec)

    spectral_lines[0].set_color("red")
            
    ax2.set_title(rf"Probability density of eigenstate $|\psi_{{{0}}}|^2$")
    ax2.set_xlabel("Position, x")
    ax2.set_ylabel("Probability density")
    ax2.set_xlim([domain[0], domain[1]])
    # plot the potential function
    potl, = ax2.plot(xgrid, V(xgrid), color="orange", zorder=5)
    # initially plot the ground state probability density on the eigenstates plot
    current_state, = ax2.plot(xgrid, energies[0] + eigenstates[:, 0] ** 2, color="blue", zorder=10)
    fill = ax2.fill_between(xgrid, energies[0], energies[0] + eigenstates[:, 0] ** 2, color="blue", alpha=0.4)
    fill_dict = {"fill": fill}

    return fig


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
    # double-well potential
    dwp = lambda x: 0.5 * (x ** 2 - 1) ** 2
    # linear ramp
    lr = lambda x, k: k * x * (x > 0)
    # morse potential
    mp = lambda x: 30 * (1 - np.exp(-1 * (x - 1))) ** 2
    # lennard-jones potential
    ljp = lambda x, eps, s: 4 * eps * ((s/x) ** 12 - (s/x) ** 6)
    # ------------------------------- #

    parser = argparse.ArgumentParser(description="Enter a one-dimensional potential")
    parser.add_argument("--potential", metavar="0.5*x**2", type=str, required=True,
                        help="Potential")
    parser.add_argument("-c", "--coeffs", nargs="+", type=float, required=False,
                        help="Superposition coefficients")
    parser.add_argument("-d", "--domain", nargs="?", required=False, help="compact interval on which the solution is found")
    parser.add_argument("-N", type=int, metavar=4, nargs="?", required=False, const=1, default=4,
                        help="Number of energy states to find")
    parser.add_argument("--densities", action='store_true')
    parser.add_argument("-v", type=bool, metavar="N", nargs="?", required=False, const=1, default=True,
                        help="Verbose mode")
    
    # get the command line arguments
    args =parser.parse_args()

    if args.coeffs:
        coeffs = np.array(args.coeffs)
        coeffs /= np.linalg.norm(coeffs)

    pot_func = lambda x: eval(args.potential)
    domain_interval = np.array(args.domain.split(","), dtype=float)

    model = Hamiltonian(potential=pot_func, domain=domain_interval)

    v, w = model.tise_solve(args.N)
    finished_plot = draw_figure(v, w, V=pot_func, domain=domain_interval)

    plt.show(block=True)