import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class StationaryEigFigure:

    def __init__(self, energies, eigenstates, V, domain, density, num_points):

        self.n = 0 # current state index
        self.num_eigenstates = len(energies)
        self.xgrid = np.linspace(domain[0], domain[1], num_points) # domain of solution
        self.V = V # potential function 
        self.spectral_lines = [] # container for holding data of spectral lines
        self.density = density

        self.energies = energies
        self.eigenstates = eigenstates
        self.current_state = None

        self.real_part = None # real part of eigenstate - only used if density option not selected
        self.im_part = None # imaginary part of eigenstate - only used if density option not selected

        self.fill = None # fill for the density - only used if the probability density option is selected 
        
        self.tot_frames = 120 # total number of frames in animation
        self.ani = None # instantiating the animation object

        # instantiate the figure
        self.fig, (self.spec_ax, self.eig_ax) = plt.subplots(nrows=1, ncols=2, figsize=(12,8),
                                                             gridspec_kw={'width_ratios': [1,4]}, sharey=True)

        self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        # self.fig.canvas.manager.set_window_title(f'Potential: {user_entered_potential}')

        # spectrum plot
        # use 1.5 as margin because the eigenstate has spatial delocalization to the extent of the well,
        # and its square must be normalized to 1. This guarantees we will never have clipping around the tops and bottoms
        top = np.max(energies) + 1.5
        self.spec_ax.set_ylabel("Energy (Hartrees)", fontsize=14)
        self.spec_ax.set_title("Spectrum", fontsize=14)
        self.spec_ax.set_ylim([min(V(self.xgrid))-1.5, top])
        # don't need the x-axis ticks and labels
        self.spec_ax.set_xticks([])
        self.spec_ax.set_xticklabels([])
        # set the initial data for the spectrum plot
        for E in self.energies:
            spec = self.spec_ax.axhline(y=E, linestyle="-", color="black")
            self.spectral_lines.append(spec)
        # change the color for the currently selected energy eigenstate
        self.spectral_lines[0].set_color("red")

        # eigenstates plot
        self.eig_ax.set_xlabel(r"Potential ($a_0$)", fontsize=14)
        self.eig_ax.set_xlim([domain[0], domain[1]])
        # plot the potential function
        self.eig_ax.plot(self.xgrid, self.V(self.xgrid), color="orange", zorder=5, label="Potential")

        # if the density option is selected, plot the probability density on the eigenstates plot (start with n=0)
        if self.density:
            self.eig_ax.set_title(rf"Probability density of eigenstate $|\psi_{{{0}}}|^2$", fontsize=14)
            self.current_state, = self.eig_ax.plot(self.xgrid, 
                                                   self.energies[0] + self.eigenstates[:, 0] ** 2, color="blue", zorder=10)
            self.fill = self.eig_ax.fill_between(self.xgrid, 
                                                 self.energies[0], self.energies[0] + self.eigenstates[:, 0] ** 2, 
                                                 color="blue", alpha=0.4)
        
        # otherwise, animate the real and imaginary parts of the eigenstate (more involved)
        else:
            self.eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{0}}}$", fontsize=14)
            self.real_part, = self.eig_ax.plot(self.xgrid, self.energies[0] + self.eigenstates[:, 0], 
                                               color="blue", zorder=10, label="Real part")
            self.im_part, = self.eig_ax.plot(self.xgrid, self.energies[0] + np.zeros_like(self.eigenstates[:, 0]), 
                                             color="pink", zorder=10, label="Imaginary part")
            self.eig_ax.legend(loc="best")
            
    
    def on_press(self, event):

        if event.key in ["left", "right"]:

            # print(f"Current state index: {n}")
            if event.key == 'left':
                inc = -1
            if event.key == 'right':
                inc = 1

            # first change the spectrum plot, meaning change the selected line's color to red 
            self.spectral_lines[self.n].set_color("black")
            # update the current state index
            self.n = (self.n + inc) % self.num_eigenstates
            self.spectral_lines[self.n].set_color("red")

            # second, change the eigenstates plot
            # plot the probability density
            if self.density:
                self.eig_ax.set_title(rf"Probability density of eigenstate $|\psi_{{{self.n}}}|^2$")
                self.current_state.set_ydata(self.energies[self.n] + self.eigenstates[:, self.n] ** 2)
                fill.remove()
                fill = self.eig_ax.fill_between(self.xgrid,
                                                self.energies[self.n], self.energies[n] + self.eigenstates[:, self.n] ** 2,
                                                color="blue", alpha=0.4)

            # animate the real and imaginary parts of the eigenstate (more involved)
            else:
                self.eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{self.n}}}$", fontsize=14)
                self.real_part.set_ydata(self.energies[self.n] + self.eigenstates[:, self.n])
                self.im_part.set_ydata(self.energies[self.n] + np.zeros_like(self.eigenstates[:, self.n]))
                self.eig_ax.legend(loc="best")

            self.fig.canvas.draw()

    
    def forward(self, frame):
        
        fundamental_freq = 1
        shift = fundamental_freq - self.energies[0]
        norm_freq = self.energies[self.n] + shift
        period = 2 * np.pi / fundamental_freq
        t = np.linspace(0, period, self.tot_frames)
        re = self.energies[self.n] + np.cos(norm_freq * t[frame]) * self.eigenstates[:, self.n]
        im = self.energies[self.n] + np.sin(norm_freq * t[frame]) * self.eigenstates[:, self.n]
        self.real_part.set_ydata(re)
        self.im_part.set_ydata(im)

        return self.real_part, self.im_part


    def begin_animate(self):
        self.ani = animation.FuncAnimation(fig=self.fig, func=self.forward, frames=self.tot_frames, interval=10, blit=True)
        plt.show()


    # method for creating and saving an animation gif that moves up the determined energies and
    # eigenstates and shows the oscillation of the real and imaginary parts
    def save_gif(self):

        return