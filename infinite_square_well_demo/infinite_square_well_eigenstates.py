import matplotlib.pyplot as plt
import numpy as np


# Constants
L = 1
m = 1
hbar = 1
n_0 = 5


def calculate_eigenstate(n, x):
    return np.sqrt(2/L)*np.sin(n*np.pi*x/L)


def calculate_allowed_energy(n):
    return (n*np.pi*hbar)**2 / (2*m*L**2)


xs = np.linspace(0, L, 10000)
ys = calculate_eigenstate(n_0, xs)
energy = calculate_allowed_energy(n_0)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')

ax.tick_params(axis='both',
               which='both',
               bottom=False,
               top=False,
               left=False,
               right=False,
               labelbottom=False,
               labelleft=False)

plt.plot(xs, ys)
plt.show()