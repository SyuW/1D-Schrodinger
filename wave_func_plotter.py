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
        
        return np.zeros(len(x))


    def solve_wavefunc(self, x, phi, dx, E):
        '''Solves for the wave function for a given eigen-energy, using Numerov's method'''

        x = x.copy()
        phi = phi.copy()

        g = 2*E*np.ones(len(x))
        f = 1 + ((dx**2)/12)*g

        for i in range(len(phi))[1:-1]:
            phi[i+1] = ((12 - 10*f[i])*phi[i] - f[i-1]*phi[i-1]) / f[i+1]
        
        return phi


    def time_evolve(self, psi, E, dt):
        '''Time evolution method for animating'''
        return


    def main(self):
        dx = 0.001
        L = 10; n = 10
        xs = np.arange(0, L, dx)

        E_0 = (n*np.pi)**2 / (2*L**2)

        phi = np.zeros(len(xs)); phi[1] = 1
        phi = self.solve_wavefunc(xs, phi, dx, E_0)
        phi /= self.normalize_state(phi, dx)
        
        plt.plot(xs, phi)
        plt.grid(True)
        plt.show()


    def __init__(self):
        self.main()


qs = QuantumSystem()