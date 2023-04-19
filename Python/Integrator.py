import numpy as np

class Integrator(object):

    """
    This temporal integrator solves an ODE of the form
    x'(t) = L*x(t)+f(x,t)
    
    Base implementation: forward Euler
    (x^{n+1}-x^{n})/dt = L*x^{n} + f(x^n,t^n)
    """

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        self._N = Nunknowns;
        self._LinearPart = LinearPart;
        self._NonLinearPart = NonlinearPart;
        self._dt = dt;
        self._iT = 0;
        alphaL = 0; alphaI = 1/dt;
        self._LinearPart.addIdentityAndFactor(alphaI,alphaL)
        
 
    def Advance(self,x0,nSteps,xprev=None): 
        """
        This is forward Euler:
        The update will be 
        x^(n+1) = x^n + dt*(L*x^n+f(x^n,t^n))
        """
        x = x0;
        for iT in range(nSteps):
            t = self._dt*self._iT;
            RHS = self._NonLinearPart.Evaluate(x,t);
            RHS += self._LinearPart.Apply(x);
            RHS += x/self._dt;
            x = RHS*self._dt;
            self._iT+=1;
        return x, xprev;

             
class BackwardEuler(Integrator):

    """
    Backward Euler:
    (x^(n+1)-x^n)/dt = L x^(n+1) + f(x^n,t^n)
    """

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        self._N = Nunknowns;
        self._LinearPart = LinearPart;
        self._NonLinearPart = NonlinearPart;
        self._dt = dt;
        self._iT = 0;
        alphaL = -1; alphaI = 1/dt;
        self._LinearPart.addIdentityAndFactor(alphaI,alphaL)
    
    def Advance(self,x0,nSteps,xprev=None): 
        """
        This is backward Euler:
        The update will be 
        x^(n+1)/dt = x^n/dt + (L*x^(n+1)+f(x^n,t^n))
        """
        x = x0;
        for iT in range(nSteps):
            t = self._dt*self._iT;
            RHS = self._NonLinearPart.Evaluate(x,t);
            RHS += x/self._dt;
            x = self._LinearPart.SolveSystem(RHS);
            self._iT+=1;
        return x, xprev;
        
class CrankNicolson(Integrator):

    """
    Crank Nicolson:
    (x^(n+1)-x^n)/dt = L (x^(n+1)/2+x^n/2) + f(3/2*x^n-1/2*x^(n-1),t^n+dt/2)
    For the nonlinear part I am combining with linear multistep in x 
    to get second order accuracy
    """

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        self._N = Nunknowns;
        self._LinearPart = LinearPart;
        self._NonLinearPart = NonlinearPart;
        self._dt = dt;
        self._iT = 0;
        alphaL = -0.5; alphaI = 1/dt;
        self._LinearPart.addIdentityAndFactor(alphaI,alphaL)
            
    def Advance(self,x0,nSteps,xprev=None): 
        """
        This is Crank Nicolson:
        The update will be 
        x^(n+1)/dt = x^n/dt + (L*x^(n+1)/2+L*x^n/2+f(x^n,t^n))
        """
        x = x0;
        if (xprev is None):
            xprev = x0;
        for iT in range(nSteps):
            t = self._dt*self._iT;
            RHS = self._NonLinearPart.Evaluate(1.5*x-0.5*xprev,t+self._dt/2);
            RHS += 0.5*self._LinearPart.Apply(x);
            RHS += x/self._dt;
            xprev = x;
            x = self._LinearPart.SolveSystem(RHS);
            self._iT+=1;
        return x, xprev;
        
