Here's my submission for the computational quantum mechanics assignment in Advanced Computational Physics (33-456). In this code, I use the shooting method for solving for the energy eigenvalues of the 1D Schroedinger equation with a harmonic oscillator potential.

We've already got an analytical solution: ( n + 0.5) * omega for n = 0, 1, 2, ...,

For this assignment, I didn't encounter any problems with the shooting method that didn't exist with the matching method as well.

# Dependencies:
Python 3, numpy
matplotlib not necessary, but you can use it for some pretty plots!

To run the code, type the following into the command line:

`python numerov.py`

You'll be prompted for some input. Below is an example:

![example](example_input.PNG)