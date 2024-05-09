import argparse
import sys

from scipy.integrate import simpson
from scipy.linalg import eigh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Hamiltonian():

    def __init__(self, potential, domain, num_points):
        """
        Constructor for the Hamiltonian class

        potential should be a callable function representing the potential energy

        """
        self.V = potential
        self.domain = domain
        self.num_points = num_points

        self.xvals = None
        self.energies = None
        self.eigenstates = None


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

        where H is the Hamiltonian defined by H = -0.5 * D^2 + V(x). 
        Note that the fact the (0,0) and (N-1, N-1) entries of the
        kinetic energy operator are equal to -2 means that Dirichlet boundary conditions are
        being applied to the endpoints of the interval.  
        """

        self.xvals = np.linspace(self.domain[0], self.domain[1], self.num_points)

        dx = self.xvals[1] - self.xvals[0]
        # off diagonals
        off = np.ones(self.num_points - 1)
        # on diagonal
        on = np.ones(self.num_points)
        # kinetic energy operator
        d_squared = 0.5 * (np.diag(off, k=1) - 2 * np.diag(on, k=0) + np.diag(off, k=-1)) / (dx ** 2)

        # Hamiltonian
        H = -d_squared + np.diag(self.V(self.xvals))

        # solve for eigenvalues/eigenstates and sort in order of increasing energy
        self.energies, self.eigenstates = eigh(a=H, subset_by_index=[0,n_eigenstates_to_find-1])

        # normalize with simpson's rule
        print("#", "-"*30, "Energy Levels", "-"*30, "#")
        for ii in range(n_eigenstates_to_find):
            self.eigenstates[:, ii] /= np.sqrt(simpson(self.eigenstates[:, ii] ** 2, x=self.xvals))
            print(f"Energy Level {ii}: {self.energies[ii]}")
        print("#", "-"*75, "#")

        return self.energies[:n_eigenstates_to_find], self.eigenstates[:, :n_eigenstates_to_find]


    # method for creating and saving an animation gif that moves up the determined energies and
    # eigenstates and shows the oscillation of the real and imaginary parts
    def make_animation(self):

        # set up the figure and axes
        fig, (spec_ax, eig_ax) = plt.subplots(nrows=1, ncols=2, figsize=(12,8),
                                              gridspec_kw={'width_ratios': [1,4]}, sharey=True)
        fig.canvas.manager.set_window_title(f"Potential: ")

        # time in-between transitions (6 secs)
        time_between = 6

        # build the gif
        frames = []

        return
    

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
                eig_ax.set_title(rf"Probability density of eigenstate $|\psi_{{{n}}}|^2$")
                current_state.set_ydata(eigenstates[:, n] ** 2 + energies[n])
                fill = fill_dict["fill"]
                fill.remove()
                fill = eig_ax.fill_between(xgrid, energies[n],
                                               energies[n] + eigenstates[:, n] ** 2,
                                               color="blue", alpha=0.4)
                fill_dict["fill"] = fill

            # animate the real and imaginary parts of the eigenstate (more involved)
            else:
                eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{n}}}$", fontsize=14)
                real_part.set_ydata(energies[n] + eigenstates[:, n])
                im_part.set_ydata(energies[n] + np.zeros_like(eigenstates[:, n]))
                eig_ax.legend(loc="best")

            fig.canvas.draw()

    # function for updating during animation loop
    def update(frame):
        n = whats_the_current_state["current"]
        # we can always shift the entire spectrum by adding a constant to the Hamiltonian
        # let us do so such that the temporal frequency of the ground state is exactly 1
        fundamental_freq = 1
        shift = fundamental_freq - energies[0]
        norm_freq = energies[n] + shift # normalized frequency
        period = 2 * np.pi / fundamental_freq
        t = np.linspace(0, period, tot_frames)
        re = energies[n] + np.cos(norm_freq * t[frame]) * eigenstates[:, n]
        im = energies[n] + np.sin(norm_freq * t[frame]) * eigenstates[:, n]
        real_part.set_ydata(re)
        im_part.set_ydata(im)

        return real_part, im_part

    # set up the figure and axes
    fig, (spec_ax, eig_ax) = plt.subplots(nrows=1, ncols=2, figsize=(12,8),
                                          gridspec_kw={'width_ratios': [1,4]}, sharey=True)
    fig.canvas.mpl_connect('key_press_event', on_press)
    fig.canvas.manager.set_window_title(f'Potential: {user_entered_potential}')

    xgrid = np.linspace(domain[0], domain[1], num_points)

    # spectrum plot
    # use 1.5 as margin because the eigenstate has spatial delocalization to the extent of the well,
    # and its square must be normalized to 1. This guarantees we will never have clipping around the tops and bottoms
    top = np.max(energies) + 1.5
    spec_ax.set_ylabel("Energy (Hartrees)", fontsize=14)
    spec_ax.set_title("Spectrum", fontsize=14)
    spec_ax.set_ylim([min(V(xgrid))-1.5, top])
    # don't need the x-axis ticks and labels
    spec_ax.set_xticks([])
    spec_ax.set_xticklabels([])

    spectral_lines = []
    for E in energies:
        spec = spec_ax.axhline(y=E, linestyle="-", color="black")
        spectral_lines.append(spec)

    # change the color for currently selected energy eigenstate
    spectral_lines[0].set_color("red")

    eig_ax.set_xlabel(r"Position ($a_0$)", fontsize=14)
    eig_ax.set_xlim([domain[0], domain[1]])
    # plot the potential function
    eig_ax.plot(xgrid, V(xgrid), color="orange", zorder=5, label="Potential")

    # plot the ground state probability density on the eigenstates plot
    if density:
        eig_ax.set_title(rf"Probability density of eigenstate $|\psi_{{{0}}}|^2$", fontsize=14)
        current_state, = eig_ax.plot(xgrid, energies[0] + eigenstates[:, 0] ** 2, color="blue", zorder=10)
        fill = eig_ax.fill_between(xgrid, energies[0], energies[0] + eigenstates[:, 0] ** 2, color="blue", alpha=0.4)
        fill_dict = {"fill": fill}

    # otherwise, animate the real and imaginary parts of the eigenstate (more involved)
    else:
        eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{0}}}$", fontsize=14)
        real_part, = eig_ax.plot(xgrid, energies[0] + eigenstates[:, 0], color="blue", zorder=10, label="Real part")
        im_part, = eig_ax.plot(xgrid, energies[0] + np.zeros_like(eigenstates[:, 0]), color="pink", zorder=10, label="Imaginary part")
        eig_ax.legend(loc="best")
        tot_frames = 120
        ani = animation.FuncAnimation(fig=fig, func=update, frames=tot_frames, interval=10, blit=True)

    plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--use_density", action='store_true', default=False,
                        help="Whether to display the probability density or the time evolving real and imaginary parts.")
    parser.add_argument("--save_gif", action='store_true', default=False,
                        help="Whether to save the animation loop as a GIF.")
    args = parser.parse_args()

    # for translating so that special functions can be parsed into corresponding numpy functions
    func_dict = {
        "sin": "np.sin",
        "cos": "np.cos",
        "sinh": "np.sinh",
        "cosh": "np.cosh",
        "tanh": "np.tanh",
        "tan": "np.tan",
        "log": "np.log",
        "log10": "np.log10",
        "ln": "np.log",
        "exp": "np.exp",
        "e": "(np.e)"
    }

    # get the user's input directly instead of through argparse and also do argument checking
    user_entered_potential = input("Enter the potential (as a function of x): ")
    
    # replace with corresponding numpy functions
    for keyword in func_dict:
        user_entered_potential = user_entered_potential.replace(keyword, func_dict[keyword])

    # now get the left and right endpoints of the domain of solution
    while True:
        try:
            left_end = input("Enter the left endpoint of the interval: ")
            left_end = float(left_end)
            break
        except ValueError:
            print("Error: left endpoint could not be converted to a real number, please try again.")
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again.")
    while True:
        try:
            right_end = input("Enter the right endpoint of the interval: ")
            right_end = float(right_end)
            break
        except ValueError:
            print("Error: right endpoint could not be converted to a real number, please try again.")
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again.")

    domain_interval = [left_end, right_end]

    # need to verify that the potential works over the domain with no singularities
    pot_func = lambda x: eval(user_entered_potential)
    try:
        pot_array = pot_func(np.linspace(*domain_interval, 100))
    except NameError:
        raise NameError("Please check the spelling/formatting in the entered potential. Exiting..")
    except Exception as e:
        print("An error occurred:", str(e), ". Please double check the spelling/formatting in the entered potential before trying again. Exiting..")
        sys.exit(1)

    # input checking for number of eigenstates to find, as well as number of grid points in discretization
    while True:
        try:
            num_eigenstates_to_find = input("Enter the amount of eigenstates you want to solve for: ")
            num_grdpts = input("Enter the number of grid points in discretization: ")
            num_eigenstates_to_find = int(num_eigenstates_to_find)
            num_grdpts = int(num_grdpts)
            if not isinstance(num_eigenstates_to_find, int):
                raise ValueError("Provided number of eigenstates to find is not an integer.")
            elif not isinstance(num_grdpts, int):
                raise ValueError("Provided number of discretization (grid-) points is not an integer.")
            elif num_eigenstates_to_find > num_grdpts:
                raise ValueError(f"The maximum number of eigenstates that can be found is the number of gridpoints: {num_grdpts}.")
            break
        except ValueError as ve:
            print(ve)
            sys.exit(1)
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again. Exiting..")
            sys.exit(1)

    announcement = f"""
    Solving for {num_eigenstates_to_find} eigenstates for the potential: {user_entered_potential} on {domain_interval} with {num_grdpts} grid points..
    """
    print(announcement)
    
    model = Hamiltonian(potential=pot_func, domain=domain_interval, num_points=num_grdpts)
    v, w = model.tise_solve(num_eigenstates_to_find)

    if args.use_density:
        print("Option '--use_density' selected. Displaying eigenstate probability densities (Born rule).")

    # create interactive animation
    animate_figure(v, w, V=pot_func, domain=domain_interval, density=args.use_density, num_points=num_grdpts)

    # give the option to user whether to save the animation or not
    if args.save_gif:
        print("Option '--save_gif' selected. Saving the animation as a GIF.")
        model.make_animation()