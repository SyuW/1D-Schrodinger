from scipy.integrate import odeint, simps, simpson
from scipy.linalg import eigh
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


class Hamiltonian():

    def __init__(self, potential, domain, num_points):
        """
        Constructor for the Hamiltonian class

        potential should be a callable function representing the potential energy

        """
        self.V = potential
        self.domain = domain
        self.num_points = num_points


    # method for solving the TISE for the Hamiltonian
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

        x_vals = np.linspace(self.domain[0], self.domain[1], self.num_points)

        dx = x_vals[1] - x_vals[0]
        off = np.ones(self.num_points - 1) # off diagonals
        on = np.ones(self.num_points) # on diagonal
        d_squared = 0.5 * (np.diag(off, k=1) - 2 * np.diag(on, k=0) + np.diag(off, k=-1)) / (dx ** 2)

        # Hamiltonian
        H = -d_squared + np.diag(self.V(x_vals))

        # solve for eigenvalues/eigenstates and sort in order of increasing energy
        energies, eigenstates = eigh(a=H, subset_by_index=[0,n_eigenstates_to_find-1])

        # normalize with simpson's rule
        for ii in range(n_eigenstates_to_find):
            eigenstates[:, ii] /= np.sqrt(simps(eigenstates[:, ii] ** 2, x=x_vals))
            print(f"Energy Level {ii}: {energies[ii]}")

        return energies[:n_eigenstates_to_find], eigenstates[:, :n_eigenstates_to_find]
    

    # method for solving the ordinary Schrodinger equation given the Hamiltonian
    def se_solve(self):
        return


def animate_figure(energies, eigenstates, V, domain, density, num_points):
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

            n = whats_the_current_state["current"]
            # print(f"Current state index: {n}")
            if event.key == 'left':
                inc = -1
            if event.key == 'right':
                inc = 1

            # first change the spectrum plot, meaning change the selected line's color to red 
            spectral_lines[n].set_color("black")
            # update the current state index
            n = (n + inc) % num_eigenstates
            whats_the_current_state["current"] = n

            spectral_lines[n].set_color("red")

            # second, change the eigenstates plot
            # plot the probability density
            if density:
                ax2.set_title(rf"Probability density of eigenstate $|\psi_{{{n}}}|^2$")
                current_state.set_ydata(eigenstates[:, n] ** 2 + energies[n])
                fill = fill_dict["fill"]
                fill.remove()
                fill = ax2.fill_between(xgrid, energies[n],
                                               energies[n] + eigenstates[:, n] ** 2,
                                               color="blue", alpha=0.4)
                fill_dict["fill"] = fill

            # animate the real and imaginary parts of the eigenstate (more involved)
            else:
                ax2.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{n}}}$")
                real_part.set_ydata(energies[n] + eigenstates[:, n])
                im_part.set_ydata(energies[n] + np.zeros_like(eigenstates[:, n]))
                ax2.legend(loc="best")

            fig.canvas.draw()

    # function for updating during animation loop
    def update(frame):
        n = whats_the_current_state["current"]
        norm_freq = energies[n] / energies[-1] # normalized frequency
        period = 2 * np.pi / (energies[0] / energies[-1])
        t = np.linspace(0, period, tot_frames)
        re = energies[n] + np.cos(norm_freq * t[frame]) * eigenstates[:, n]
        im = energies[n] + np.sin(norm_freq * t[frame]) * eigenstates[:, n]
        real_part.set_ydata(re)
        im_part.set_ydata(im)

        return real_part, im_part

    # set up the figure and axes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,8), 
                                   gridspec_kw={'width_ratios': [1,4]}, sharey=True)
    fig.canvas.mpl_connect('key_press_event', on_press)
    fig.canvas.manager.set_window_title(f'Potential: {args.potential}')

    xgrid = np.linspace(domain[0], domain[1], num_points)

    # spectrum plot
    # use 1.5 because the eigenstate has spatial delocalization to the extent of the well, and must be normalized to 1
    top = np.max(energies) + 1.5
    ax1.set_ylabel("Energy, E")
    ax1.set_title("Spectrum")
    ax1.set_ylim([min(V(xgrid))-1.5, top])
    # don't need the x-axis ticks and labels
    ax1.set_xticks([])
    ax1.set_xticklabels([])

    spectral_lines = []
    for E in energies:
        spec = ax1.axhline(y=E, linestyle="--", color="black")
        spectral_lines.append(spec)

    spectral_lines[0].set_color("red")

    ax2.set_xlabel("Position, x")
    ax2.set_xlim([domain[0], domain[1]])
    # plot the potential function
    ax2.plot(xgrid, V(xgrid), color="orange", zorder=5)

    # plot the ground state probability density on the eigenstates plot
    if density:
        ax2.set_title(rf"Probability density of eigenstate $|\psi_{{{0}}}|^2$")
        ax2.set_ylabel("Probability density")
        current_state, = ax2.plot(xgrid, energies[0] + eigenstates[:, 0] ** 2, color="blue", zorder=10)
        fill = ax2.fill_between(xgrid, energies[0], energies[0] + eigenstates[:, 0] ** 2, color="blue", alpha=0.4)
        fill_dict = {"fill": fill}

    # otherwise, animate the real and imaginary parts of the eigenstate (more involved)
    else:
        ax2.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{0}}}$")
        real_part, = ax2.plot(xgrid, energies[0] + eigenstates[:, 0], color="blue", zorder=10, label="Real part")
        im_part, = ax2.plot(xgrid, energies[0] + np.zeros_like(eigenstates[:, 0]), color="pink", zorder=10, label="Imaginary part")
        ax2.legend(loc="best")
        tot_frames = 120
        ani = animation.FuncAnimation(fig=fig, func=update, frames=tot_frames, interval=10, blit=True)

    plt.show()


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
    dwp = lambda x, a: 0.5 * (x ** 2 - a ** 2) ** 2
    # linear ramp
    lr = lambda x, k: k * x * (x > 0)
    # morse potential
    mp = lambda x: 30 * (1 - np.exp(-1 * (x - 1))) ** 2
    # lennard-jones potential
    ljp = lambda x, eps, s: 4 * eps * ((s/x) ** 12 - (s/x) ** 6)
    # asymmetric double well potential
    adwp = lambda x: 0.5 * ((x - 2) ** 2) * ((x + 2) ** 2 + 1) 
    # ------------------------------- #

    parser = argparse.ArgumentParser(description="Enter a one-dimensional potential")
    parser.add_argument("--potential", metavar="0.5*x**2", type=str, required=True,
                        help="Potential")
    parser.add_argument("-c", "--coeffs", nargs="+", type=float, required=False,
                        help="Superposition coefficients")
    parser.add_argument("-d", "--domain", nargs="?", required=False, help="compact interval on which the solution is found")
    parser.add_argument("-N", type=int, metavar=4, nargs="?", required=False, const=1, default=4,
                        help="Number of energy states to find")
    parser.add_argument("--num_grdpts", type=int, nargs="?", required=False, default=1000, help="Number of grid points to use")
    parser.add_argument("--use_density", action='store_true', default=False)
    parser.add_argument("-v", action='store_true', required=False, default=True,
                        help="Verbose mode")
    
    # get the command line arguments
    args =parser.parse_args()

    if args.coeffs:
        coeffs = np.array(args.coeffs)
        coeffs /= np.linalg.norm(coeffs)

    pot_func = lambda x: eval(args.potential)
    domain_interval = np.array(args.domain.split(","), dtype=float)

    model = Hamiltonian(potential=pot_func, domain=domain_interval, num_points=args.num_grdpts)

    v, w = model.tise_solve(args.N)
    animate_figure(v, w, V=pot_func, domain=domain_interval, density=args.use_density, num_points=args.num_grdpts)