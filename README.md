# 1D-Schrodinger
A script for solving the 1D time-independent Schrodinger equation for arbitrary user-specified potentials, giving the quantized energy levels of a trapped quantum particle.

Produces an ✨interative✨ animated plot for viewing the energy spectrum and eigenstates. Simply use the left and right arrow keys to toggle between the different excited modes of the particle.  

# Usage
Simply run using the command: ```python schrodinger.py```

You will be prompted to enter the potential, left endpoint, right endpoint, number of eigenstates to solve for, and number of gridpoints to use.

An example is shown below for an asymmetric double well potential of the form ```0.5 * ((x - 2) ** 2) * ((x + 2) ** 2 + 1)```. Note the small energy gap between the 3rd and 4th energy levels (contributions from quantum tunnelling AKA instantons breaks the degeneracy):

![Asymmetric double well](https://github.com/SyuW/1D-Schrodinger/blob/master/demos/double_well.gif)

## Useful Links
A collection of resources I referenced and took inspiration from throughout this work. Check them out!
- http://www.pas.rochester.edu/~tobin/notebook/2009/02/12/html/eigfunc.html. Uses the diagonalization approach for solving the 1D TISE. Finds the eigenvalues/eigenstates iteratively by applying the Hamiltonian as an iterated map acting on an initial random vector. 
- https://github.com/FelixDesrochers/Numerov. Has an implementation of Numerov's method for solving the 1D TISE for bound states as well as pretty figures.
- Professor Qijing Zheng's Website: http://staff.ustc.edu.cn/~zqj/post/. Has a couple of useful posts/code for solving the 1D SE and TISE, also has nice writing about other topics in computational condensed matter/solid state physics if you're interested in that.
- http://physics.unipune.ac.in/~phyed/23.1/23.1_computation.pdf. Goes into some numerical instability issues with the Numerov method, such as when integrating too far into the classically forbidden region of a potential. Finds a forwards and backwards integrated solution and then imposes matching condition based off the continuity of second derivative.
- https://jakevdp.github.io/blog/2012/09/05/quantum-python/. Has a tutorial for using a split step Fourier method for solving the SE.
