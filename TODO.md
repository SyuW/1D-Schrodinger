## Features to add
- GIF creation
- Don't hardcode the upper and lower padding as 1.5: use classical turning points of highest energy level + normalization of probability to set bounded area
- Construct time-varying superposition of states: would be interesting to see quantum tunnelling through central barrier in the double well potential
- Evolve a quantum state in time starting from an initial state
- Extend to time-varying potentials: visualize the effects of geometric phase 
- Embed the solver into a webpage: need a port to Javascript & WebGL for this??

## More distant future goals
- Add support for 2D and 3D potentials
- Extend solver to relativistic wave equations such Dirac or Klein-Gordon
- Extend solver to multiple particles: examine the effects of quantum statistics
- Extend solver to nonlinear variants of Schrodinger equation, such as the Gross-Pitaevskii equation used for studying trapped ultracold atoms
- Visualization of states

## Refactoring
- Wrap figure creation and animation into a class instead of having nested functions

## Bugs
- Figure out how to make animation of excited states into a perfect loop while still depicting the oscillation speed 
- For the double well potential, sometimes the ground state does not have parity symmetry
