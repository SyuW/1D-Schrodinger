# 1D-Schrodinger
Experimenting with methods for solving the Schrodinger equation in one dimension, both time-dependent and time-independent variants.
The motivation is to gain intuition for the dynamics of bound quantum states in 1D and also pick up cool tricks along the way for
solving differential equations numerically 😃!

## First phase: Solve the 1D TISE for arbitrary potential ##
One way to construct a time-varying quantum state is through superpositions of eigenfunctions of the system's Hamiltonian operator, which form a complete orthogonal basis. Hence,
we have to solve the 1D time independent Schrodinger equation (TISE), a second order ODE.

Some methods I've tried:
- Numerov's algorithm
- Converting to a system of 1st order ODEs
- Diagonalization approach

Of the methods I've tried, **diagonalization** has been the most numerically stable as well as being pretty straightforward to implement.
This is the approach I'll mostly likely be using moving forward for solving the 1D TISE.  

## Second phase: Time evolution ##
To be continued.

## Third phase ##
To be continued.

# Useful Links
A collection of resources I referenced and took inspiration from throughout this work. Check them out!
Professor Qijing Zheng's Website: http://staff.ustc.edu.cn/~zqj/post/. Has some useful
