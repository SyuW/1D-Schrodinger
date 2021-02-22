import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


class QuantumSystem:
    '''
    1D Quantum Mechanical System Class
    '''

    def normalize_state(self, data):
        '''Calculate the normalization factor of a quantum state'''

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


    def solve_TISE(self, dx):
        '''Solves the TISE for a given eigen-energy, using Numerov's method'''

        x = self.x.copy()
        V = self.calculate_potential(x)

        tolerance = 10**(-7)
        boundary_value = 0

        def integrator(energy, pot):
            phi = self.phi.copy()
            g = (2*self.m/self.hbar**2)*(energy - pot)
            f = 1 + ((dx**2)/12)*g
            for i in range(len(phi))[1:-1]:
                phi[i+1] = ((12 - 10*f[i])*phi[i] - f[i-1]*phi[i-1]) / f[i+1]
            return phi

        def count_nodes_of_func(f):
            node_count = 0
            for i, val in enumerate(f[:-2]):
                if f[i+1] * f[i+2] <= 0: # Change of sign indicates presence of node
                    node_count += 1
            return node_count

        def bisection_method(nodes_wanted, E_1, E_2, V=V):

            E_trial = (E_2 + E_1) / 2 # Set energy guess to average initially
            test = integrator(energy=E_trial, pot=V)
            node_count = count_nodes_of_func(f=test) # The number of nodes + 1 gives corresponding n-th energy eigenvalue
            while node_count != nodes_wanted:
                if node_count > nodes_wanted: # The energy guess was too big (overshot)
                    E_2 = (E_1 + E_2) / 2
                elif node_count < nodes_wanted: # The energy guess was too small (undershot)
                    E_1 = (E_1 + E_2) / 2
                E_trial = (E_1 + E_2) / 2
                test = integrator(energy=E_trial, pot=V)
                node_count = count_nodes_of_func(f=test)

            return E_trial

        def construct_energy_guess(energy_number, E_lower, E_upper, V):

            modified_upper = bisection_method(nodes_wanted=energy_number, E_1=E_lower, E_2=E_upper, V=V)
            modified_lower = bisection_method(nodes_wanted=energy_number-1, E_1=E_lower, E_2=modified_upper, V=V)

            return modified_lower, modified_upper
            
        E_lower_limit = 0; E_upper_limit = 2000
        energy_no = 19

        lower, upper = construct_energy_guess(energy_no, E_lower_limit, E_upper_limit, V)
        
        E_trial = (upper + lower) / 2
        print(f"Lower limit: {lower}"); print(f"Upper limit: {upper}")
        phi = integrator(energy=E_trial, pot=V)

        while abs(phi[-1] - boundary_value) > tolerance:
            if energy_no % 2 != 0: # Odd energy number will give an even number of nodes
                #print(f"{E_trial}")
                if (phi[-1] - boundary_value) < 0: # If negative at end, the energy eigenvalue guess is too high
                    #print(f"Too high: {phi[-1]}")
                    upper = (upper + lower) / 2
                elif (phi[-1] - boundary_value) > 0: # If positive at end, the energy eigenvalue guess is too low
                    #print(f"Too low: {phi[-1]}")
                    lower = (upper + lower) / 2

            else: # Even energy number gives odd number of nodes
                if (phi[-1] - boundary_value) > 0: # If positive at end, the energy eigenvalue guess is too high
                    #print(f"Too high: {phi[-1]}")
                    upper = (upper + lower) / 2
                elif (phi[-1] - boundary_value) < 0: # If negative at end, the energy eigenvalue guess is too low
                    #print(f"Too low: {phi[-1]}")
                    lower = (upper + lower) / 2

            E_trial = (upper + lower) / 2
            phi = integrator(energy=E_trial, pot=V)

        print(f"Found {energy_no}th energy eigenvalue: {E_trial}")

        return phi
        


    def __init__(self):
        
        L = 1; n = 1
        self.dx = 0.0001

        self.hbar = 1
        self.m = 1
        self.ang_freq = 1

        self.E_n = (n*np.pi)**2 / (2*L**2) #infinite well

        self.x = np.arange(0, L, self.dx)
        self.time_elapsed = 0

        self.phi = np.zeros(len(self.x))
        self.phi[1] = 0.00001

        self.phi = self.solve_TISE(self.dx)
        self.phi /= self.normalize_state(self.phi)


def main():
    qs = QuantumSystem()

    y_bound = max(qs.phi) + 1
    fig = plt.figure()
    ax  = plt.axes(xlim=(min(qs.x)-1, max(qs.x)+1), ylim=(-y_bound, y_bound))
    ax.grid()
    line, = ax.plot(qs.x, qs.phi)
    plt.show()


main()