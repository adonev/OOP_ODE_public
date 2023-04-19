import numpy as np
from Integrator import Integrator, BackwardEuler, CrankNicolson
from NonStiffIntegrand import NonStiffIntegrand
from LinearOperator import DenseLinearOperator, SparseLinearOperator
from math import sin, cos

"""
Solve the ODE 
dx/dt = A*x + B, where
B(x,t) is given by NonStiffPart below
"""
def NonStiffPart(x,t):
    Part1 = 0.2e0*x[0]*x[1]+0.125142e1*cos(t)-0.911111e0*sin(t)-0.2e0*sin(t)*cos(t);
    Part2 = 0.3e0*x[0]**2-0.1e0*x[1]**2-0.748585e0*sin(t)-0.288888e0*cos(t)-0.3e0*sin(t)**2+0.1e0*cos(t)**2;
    return np.reshape(np.array([Part1,Part2]),(2,1));

eig = 10;
A = -np.mat([[0.222222e-1+0.888889e0*eig,0.628539e-1-0.314269e0*eig],[0.628536e-1-0.314269e0*eig,0.177777e0+0.111111e0*eig]])
N = 2;
dt = 0.004
tf = 2;
x0=np.ones((N,1))
LinearPart = DenseLinearOperator(N);
LinearPart.setOperator(A,invertible=True);
ExplicPart = NonStiffIntegrand(NonStiffPart);

# Code to illustrate how to call ODE solvers
TempInt = CrankNicolson(N,LinearPart,ExplicPart,dt);
#TempInt = BackwardEuler(N,LinearPart,ExplicPart,dt);
#TempInt = Integrator(N,LinearPart,ExplicPart,dt);
nSteps = int(tf/dt+1e-6);
xFinal, _ =TempInt.Advance(x0,nSteps);
xTrue=np.array([[-0.253904123470440],[-0.469937407560826]]) # From matlab (accurate to 10 digits)
print(np.linalg.norm(xFinal-xTrue))

