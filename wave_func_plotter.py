import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


class QuantumSystem:
    '''
    1D Quantum Mechanical System Class
    '''

    def normalize_state(self, data):
        '''Calculate the normalization factor of quantum state'''

        data = (data.copy())**2

        # Use Simpson's Rule
        running_sum = data[0] + data[-1]
        for i in range(len(data)):
            if i % 2 == 0:
                running_sum += 4*data[i]
            else:
                running_sum += 2*data[i]

        return np.sqrt((self.dx/3)*running_sum)


    def calculate_potential(self, x):
        '''Potential Function'''

        x = x.copy()
        V = np.zeros(len(x))

        for i, e in enumerate(x):
            if abs(e) >= 3:
                V[i] = 0
            else:
                continue
        
        return V

    
    def determine_energies(self):
        N_x = 1000
        x = np.linspace(0, 10, N_x)

        middle_diag = np.full(N_x, -2)
        up_down_diag = np.full(N_x-1, 1)

        V_ij = np.diag(self.calculate_potential(x))
        K_ij = (np.diag(up_down_diag, 1) + np.diag(up_down_diag, -1) + np.diag(middle_diag)) # / self.dx**2

        H_ij = -K_ij + V_ij

        print(H_ij); print(np.size(H_ij))

        return np.linalg.eig(H_ij)


    def solve_TISE(self, dx):
        '''Solves the TISE for a given eigen-energy, using Numerov's method'''

        x = self.x.copy()
        phi = self.phi.copy()

        V = self.calculate_potential(x)

        g = (2*self.m/self.hbar**2)*(self.E_n - V)
        f = 1 + ((dx**2)/12)*g

        for i in range(len(phi))[1:-1]:
            phi[i+1] = ((12 - 10*f[i])*phi[i] - f[i-1]*phi[i-1]) / f[i+1]
        
        return phi


    def __init__(self):
        
        L = 10; n = 3
        self.dx = 0.001

        self.hbar = 1
        self.m = 1
        self.ang_freq = 2

        self.E_n = (n*np.pi)**2 / (2*L**2) #infinite well

        self.x = np.arange(0, L, self.dx)
        self.time_elapsed = 0

        self.phi = np.zeros(len(self.x))
        self.phi[1] = 1

        self.phi = self.solve_TISE(self.dx)
        self.phi /= self.normalize_state(self.phi)


def main():
    qs = QuantumSystem()

    y_bound = max(qs.phi) + 1
    fig = plt.figure()
    ax  = plt.axes(xlim=(min(qs.x)-1, max(qs.x)+1), ylim=(-y_bound, y_bound))
    ax.grid()
    line, = ax.plot([], [], lw=2)

    dt = (1./30)


    def init():
        line.set_data([], [])
        return line,


    def animate(i):
        x = qs.x
        y = qs.phi * np.cos(qs.ang_freq * dt * i)
        line.set_data(x, y)
        return line,


    from time import time
    t0 = time()
    animate(0)
    t1 = time()
    interval = 1000 * dt - (t1 - t0)

    frame_no = int((2*np.pi/qs.ang_freq) * (1 / dt)) + 1

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=frame_no, interval=interval, blit=True)

    plt.show()


def test():
    qs = QuantumSystem()
    w, v = qs.determine_energies()

    print(w)
    print(np.linalg.eig(np.diag(np.arange(3))))


test()