import numpy as np
from scipy.integrate import simpson
from scipy.linalg import eigh


class SingleParticleHamiltonian:

    def __init__(self, V, domain, num_points):
        """
        Constructor for the Hamiltonian class

        potential should be a callable function representing the potential energy

        """
        self.V = V
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
    

    # method for solving the ordinary Schrodinger equation given the Hamiltonian
    def se_solve(self):
        return


class TwoParticleHamiltonian:

    def __init__(self) -> None:
        pass
        

class ThreeParticleHamiltonian:

    def __init__(self) -> None:
        pass