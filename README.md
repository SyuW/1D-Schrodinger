# 1D-Schrodinger
Experimenting with methods for solving the Schrodinger equation in one dimension, both time-dependent and time-independent variants.
The motivation is to gain intuition for the dynamics of bound quantum states in 1D and also to pick up cool tricks along the way for
solving differential equations numerically 😃!

## First phase: Solve the 1D TISE for arbitrary potential ##
One way to construct a time-varying quantum state is through superpositions of eigenstates of the system's Hamiltonian operator, which form a complete orthogonal basis. Looking qualitatively at how the energy eigenstate probability densities depend on potential shape is also pretty instructive. Hence,
we need to solve the 1D time independent Schrodinger equation (TISE), a second order ODE.

Some methods I've tried:
- Numerov's algorithm
- Converting to a system of 1st order ODEs
- Diagonalization approach

Of the methods I've tried, **diagonalization** has been the most numerically stable as well as being pretty straightforward to implement.
This is the approach I'll mostly likely be using moving forward for solving the 1D TISE.

![Example 1](samples/linear_potential_well.png?raw=true "Linear Potential Well")

## Second phase: Time evolution ##
If the energy levels considered are smaller than the potential well depth, they are discrete, allowing us to write the time-varying state as a sum over energy eigenstates:

$$
  \psi(x,t) = \sum_{n}c_n(0)e^{iE_n{t}/\hbar}\phi_n(x)
$$

`schrodinger.py` allows you to provide coefficients $c_n(0)$. 

## Third phase: GUI? ##
In the future I intend to develop a GUI for the program instead of requiring command line input. I'm still trying to work out the layout of the interface.

# Useful Links
A collection of resources I referenced and took inspiration from throughout this work. Check them out!
- http://www.pas.rochester.edu/~tobin/notebook/2009/02/12/html/eigfunc.html. Uses the diagonalization approach for solving the 1D TISE. Finds the eigenvalues/eigenstates iteratively by applying the Hamiltonian as an iterated map acting on a random vector. 
- https://github.com/FelixDesrochers/Numerov. Has an implementation of Numerov's method for solving the 1D TISE for bound states as well as pretty figures.
- Professor Qijing Zheng's Website: http://staff.ustc.edu.cn/~zqj/post/. Has a couple of useful posts/code for solving the 1D SE, also has nice writing about other topics in computational condensed matter/solid state physics if you're interested in that.
- http://physics.unipune.ac.in/~phyed/23.1/23.1_computation.pdf. Goes into some numerical instability issues with the Numerov method, such as when integrating too far into the classically forbidden region of a potential. Finds a forwards and backwards integrated solution and then imposes matching condition based off the continuity of second derivative.
- https://jakevdp.github.io/blog/2012/09/05/quantum-python/. Has a tutorial for using a split step Fourier method for solving the SE.
