# 1D-Schrodinger
Experimenting with methods for solving the Schrodinger equation in one dimension, both time-dependent and time-independent variants.
The motivation is to gain intuition for the dynamics of bound quantum states in 1D and also to pick up cool tricks along the way for
solving differential equations numerically 😃!

## First phase: Solve the 1D TISE for arbitrary potential ##
One way to construct a time-varying quantum state is through superpositions of eigenstates of the system's Hamiltonian operator, which form a complete orthogonal basis. Looking qualitatively at how the energy eigenstates depend on potential shape is also pretty instructive. Hence,
we need to solve the 1D time independent Schrodinger equation (TISE), a second order ODE.

Some methods I've tried:
- Numerov's algorithm
- Converting to a system of 1st order ODEs
- Diagonalization approach

Of the methods I've tried, **diagonalization** has been the most numerically stable as well as being pretty straightforward to implement.
This is the approach I'll mostly likely be using moving forward for solving the 1D TISE.  

## Second phase: Time evolution ##
To be continued.

## Third phase: GUI? ##
To be continued.

# Useful Links
A collection of resources I referenced and took inspiration from throughout this work. Check them out!
- Numerical quantum mechanics by Tobin Fricke: http://www.pas.rochester.edu/~tobin/notebook/2009/02/12/html/eigfunc.html. Uses the diagonalization approach for solving the 1D TISE. Finds the eigenvalues/eigenstates iteratively by applying the Hamiltonian as an iterated map acting on a random vector. 
- https://github.com/FelixDesrochers/Numerov. Has an implementation of Numerov's method for solving the 1D TISE for bound states as well as pretty figures.
- Professor Qijing Zheng's Website: http://staff.ustc.edu.cn/~zqj/post/. Has a couple of useful posts/code for solving the 1D SE, also has nice writing about other topics in computational condensed matter/solid state physics if you're interested in that.
- Bound State of One Dimensional Potential by Numerov Method by D G Kenhere: http://physics.unipune.ac.in/~phyed/23.1/23.1_computation.pdf. Goes into some numerical instability issues with the Numerov method, such as when integrating too far into the classically forbidden region of a potential.
