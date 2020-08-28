import numpy as np
import matplotlib.pyplot as plt


def normalize_state(data, dx):
    data = (data.copy())**2

    running_sum = data[0] + data[-1]
    for i in range(len(data)):
        if i % 2 == 0:
            running_sum += 4*data[i]
        else:
            running_sum += 2*data[i]

    return np.sqrt((dx/3)*running_sum)
    

def calculate_potential(x):
    x = x.copy()
    
    return np.zeros(len(x))


def numerov_method(x, phi, dx, E):
    x = x.copy()
    phi = phi.copy()

    g = 2*E*np.ones(len(x))
    f = 1 + ((dx**2)/12)*g

    for i in np.arange(len(phi))[1:-1]:
        phi[i+1] = ((12 - 10*f[i])*phi[i] - f[i-1]*phi[i-1]) / f[i+1]
    
    return phi


def main():
    dx = 0.001
    L = 10
    n = 10
    xs = np.arange(0, L, dx)
    E_0 = (n*np.pi)**2 / (2*L**2)

    phi = np.zeros(len(xs)); phi[1] = 1
    phi = numerov_method(xs, phi, dx, E_0)
    phi /= normalize_state(phi, dx)
    
    plt.plot(xs, phi)
    plt.grid(True)
    plt.show()


main()