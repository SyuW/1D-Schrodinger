import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


class QuantumSystem:
    '''
    1D Quantum Mechanical System Class
    '''

    def normalize_state(self, data, dx):
        '''Calculate the normalization factor of quantum state'''

        data = (data.copy())**2

        # Use Simpson's Rule
        running_sum = data[0] + data[-1]
        for i in range(len(data)):
            if i % 2 == 0:
                running_sum += 4*data[i]
            else:
                running_sum += 2*data[i]

        return np.sqrt((dx/3)*running_sum)


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


    def solve_TISE(self, x, phi, dx, E):
        '''Solves the TISE for a given eigen-energy, using Numerov's method'''

        x = x.copy()
        phi = phi.copy()

        V = self.calculate_potential(x)

        g = (2*self.m/self.hbar**2)*(E - V)
        f = 1 + ((dx**2)/12)*g

        for i in range(len(phi))[1:-1]:
            phi[i+1] = ((12 - 10*f[i])*phi[i] - f[i-1]*phi[i-1]) / f[i+1]
        
        return phi


    def __init__(self):
        dx = 0.001
        L = 10; n = 3

        self.hbar = 1
        self.m = 1
        self.ang_freq = 2

        self.E_n = (n*np.pi)**2 / (2*L**2) #infinite well

        self.x = np.arange(0, L, dx)
        self.time_elapsed = 0

        self.phi = np.zeros(len(self.x))
        self.phi[1] = 1

        self.phi = self.solve_wavefunc(self.x, self.phi, dx, self.E_n)
        self.phi /= self.normalize_state(self.phi, dx)


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