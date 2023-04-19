import numpy as np

class NonStiffIntegrand(object):
# This class can have many other methods like
# Solve (df/dx+alpha*I)x=b
# or at least to return df/dx(x) as a linear operator
# if we want to make this a more general "ODEIntegrand" rather than specifically "NonStiff" integrand
# See the Matlab codes for an example of this
    
    def __init__(self,function):
        self._function = function;
    
    def Evaluate(self,x,t):
        return self._function(x,t);
        
