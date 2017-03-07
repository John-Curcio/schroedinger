Here's my submission for a project in my Advanced Computational Physics class (33-456). In this code, I use the Numerov method, together with the shooting method, for solving the energy eigenvalues of the following 1D Schroedinger equation:

![equation](equation.jpeg)

Omega is an input parameter proportional to the frequency.

The equation above describes a simple harmonic oscillator potential, for which we've already got an analytical solution. But the numerical method here gets the same result, and will work for an arbitrary external potential (so long as it's "well-behaved").
