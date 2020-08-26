import numpy as np
import matplotlib.pyplot as plt


def determine_eigenvalues(xs):
    pass


def calculate_potential(x):
    # Potential function
    if abs(x) >= 1:
        value = 100
    else:
        value = 0
    return value


def calculate_curvature(phi_x, V):
    return 2*coeff*(V - E)*phi_x


def main():
    dx = 0.001
    x_points = np.arange(-1.5, 1.5, dx)

    for x in xs:
        # First determine curvature at current x
        curva_x    = calculate_curvature(eigenstate[-1], calculate_potential(x))
        curva_x_dx = calculate_curvature(eigenstate[-1], calculate_potential(x + dx))

        # Update eigenstate at x + dx
        eigenstate += [eigenstate[-1] + eigenstate_slope[-1]*dx + 0.5*curva_x*(dx**2)] 

        # Update eigenstate slope at x + dx
        eigenstate_slope += [eigenstate_slope[-1] + 0.5*(curva_x + curva_x_dx)*dx]
    
    plt.plot(xs, np.array(eigenstate[:-1]))
    plt.show()