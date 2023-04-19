# OOP_ODE_public
Example suite of ODE solvers with object-oriented design.
We provide codes in Matlab, Python, and Fortran (not yet finished), see the corresponding subdirectories
The codes show how to design flexible, extensible, and reusable implicit-explicit (or exponential-explicit) ODE solvers.
These codes are provided as examples only. They are neither complete nor efficient as written.

The Matlab codes were written by Mariya Savinov, the python ones by Ondrej Maxian, and Fortran ones by Aleksandar Donev,
all at the Courant Institute of Mathematics at NYU at the time of writing.

The ODEs are assumed to have the form:

dx/dt = A*x + B(x,t)

where A is a LinearOperator and B is a nonlinear Integrand. 

In the default implementation of an Integrator, we use forward Euler for the nonlinear term and backward Euler for the linear term.
Also implemented in this public release are several IMEX one and two step RK schemes. 

We test the codes out on an ODE system with two variables, constructed to have the solution [sin(t),cos(t)], see the Maple subdirectory.
The parameter max_eig controls is the L2 conditioning number of A. The nonlinear term is not stiff.
We also illustrate how to solve a nonlinear diffusion equation in 1d with homogeneous DirichletBCs.

u\_t = d*u\_xx + u^2

which is a toy model of ignition. The solution can blow up in finite time if the initial data is large.

