import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class StationaryEigFigure:

    def __init__(self, energies, eigenstates, Vname, domain, density, num_points, darkmode):

        self.n = 0 # current state index
        self.num_eigenstates = len(energies)
        self.xgrid = np.linspace(domain[0], domain[1], num_points) # domain of solution
        self.V = lambda x: eval(Vname) # potential function 
        self.spectral_lines = [] # container for holding data of spectral lines
        self.density = density

        self.energies = energies
        self.eigenstates = eigenstates
        self.current_state = None
        self.fundamental_freq = 1

        # only used if density option is not selected: dynamic animation
        self.real_part = None # real part of eigenstate
        self.im_part = None # imaginary part of eigenstate
        self.tot_frames = 120 # total number of frames in animation
        self.ani = None # instantiating the animation object
        
        self.t = np.linspace(0, 2 * np.pi / self.fundamental_freq, self.tot_frames)

        # only used if the probability density option is selected
        self.fill = None # fill for the density 

        # toggling dark mode
        if darkmode:
            plt.style.use("dark_background")
            self.normal_spec_color = "white"
        else:
            self.normal_spec_color = "black"

        # instantiate the figure
        self.fig, (self.spec_ax, self.eig_ax) = plt.subplots(nrows=1, ncols=2, figsize=(12,8),
                                                             gridspec_kw={'width_ratios': [1,4]}, sharey=True)

        self.fig.canvas.mpl_connect('key_press_event', self.on_press)


        self.fig.canvas.manager.set_window_title(f'Potential: {Vname}')

        # spectrum plot
        # use 1.5 as margin because the eigenstate has spatial delocalization to the extent of the well,
        # and its square must be normalized to 1. This guarantees we will never have clipping around the tops and bottoms
        top = np.max(energies) + 1.5
        self.spec_ax.set_ylabel("Energy (Hartrees)", fontsize=14)
        self.spec_ax.set_title("Spectrum", fontsize=14)
        self.spec_ax.set_ylim([min(self.V(self.xgrid))-1.5, top])
        # don't need the x-axis ticks and labels
        self.spec_ax.set_xticks([])
        self.spec_ax.set_xticklabels([])
        # set the initial data for the spectrum plot
        for E in self.energies:
            spec = self.spec_ax.axhline(y=E, linestyle="-", color=self.normal_spec_color)
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

            # first change the spectrum plot, which means to reset the current spectral line
            # to its default color and switch the next spectral line to red
            self.spectral_lines[self.n].set_color(self.normal_spec_color)
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
                                                self.energies[self.n],
                                                self.energies[self.n] + self.eigenstates[:, self.n] ** 2,
                                                color="blue", alpha=0.4)

            # animate the real and imaginary parts of the eigenstate (more involved)
            else:
                self.eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{self.n}}}$", fontsize=14)
                self.real_part.set_ydata(self.energies[self.n] + self.eigenstates[:, self.n])
                self.im_part.set_ydata(self.energies[self.n] + np.zeros_like(self.eigenstates[:, self.n]))
                self.eig_ax.legend(loc="best")

            self.fig.canvas.draw()

    
    def forward(self, frame):
        
        shift = self.fundamental_freq - self.energies[0]
        norm_freq = self.energies[self.n] + shift
        re = self.energies[self.n] + np.cos(norm_freq * self.t[frame]) * self.eigenstates[:, self.n]
        im = self.energies[self.n] + np.sin(norm_freq * self.t[frame]) * self.eigenstates[:, self.n]
        self.real_part.set_ydata(re)
        self.im_part.set_ydata(im)

        return self.real_part, self.im_part


    def begin_interactive_animate(self):
        ani = animation.FuncAnimation(fig=self.fig, func=self.forward, frames=self.tot_frames, interval=10, blit=True)
        plt.show()


    # create an animation suitable for saving to a GIF
    def gif_animate(self, frame):

        frames_per_level = self.seconds_per_eigenstate * self.fps
        # frames_per_level = 120
        n = frame // frames_per_level # state index

        shift = self.fundamental_freq - self.energies[0]
        norm_freq = self.energies[n] + shift
        t = np.linspace(0, 2 * np.pi / self.fundamental_freq, frames_per_level) # array of time points in oscillation

        # change the spectrum plot
        self.spectral_lines[n].set_color("red")
        for i, line in enumerate(self.spectral_lines):
            if i != n:
                line.set_color(self.normal_spec_color)
        
        # change the eigenstates plot 
        self.eig_ax.set_title(rf"Real and imaginary parts of eigenstate $\psi_{{{n}}}$", fontsize=14)
        re = self.energies[n] + self.eigenstates[:, n] * np.cos(norm_freq * t[frame % frames_per_level])
        im = self.energies[n] + self.eigenstates[:, n] * np.sin(norm_freq * t[frame % frames_per_level])
        self.real_part.set_ydata(re)
        self.im_part.set_ydata(im)
        self.eig_ax.legend(loc="best")

        return self.real_part, self.im_part


    # method for creating and saving an animation gif that moves up the determined energies and
    # eigenstates and shows the oscillation of the real and imaginary parts
    def make_gif(self, max_eigenstates=None, fps=30, seconds_per_eigenstate=2):

        self.fps = fps
        self.seconds_per_eigenstate = seconds_per_eigenstate

        if not max_eigenstates or (self.num_eigenstates < max_eigenstates):
            max_eigenstates = self.num_eigenstates

        gif_total_frames = fps * seconds_per_eigenstate * max_eigenstates

        ani = animation.FuncAnimation(fig=self.fig, func=self.gif_animate, 
                                               frames=gif_total_frames, 
                                               interval=1/fps * 1000, blit=True)

        ani.save(filename="../gifs/animation.gif", writer="Pillow")

        plt.show()