import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import eye as speye
from scipy.linalg import cho_factor, lu_factor, cho_solve, lu_solve, expm
from scipy.sparse.linalg import splu as splu
from scipy.sparse.linalg import expm as sp_expm

"""
This class describes a linear operator L which acts on a vector x
"""
class LinearOperator(object): 

    """
    Abstract class is just scalars. 
    """
    def __init__(self,Nunknown):
        """
        Initialize with the number of unknowns. Initialize memory for the matrix representation of L
        """
        self._N = Nunknown;
        self._L = 0;
    
    def setOperator(self,a,invertible=True):
        """
        Set the linear operator A. 
        In this abstract class, this is just a scalar
        Invertible = boolean to say whether a is invertible or not
        """
        self._L = a;
        self._Invertible = invertible;
        
    def PreCalcExponential(self,dt):
        """
        Precompute the exponential e^(A dt), where dt is 
        the time step size
        """
        self._ExpL = np.exp(self._L*dt);
        
    def ApplyExponential(self,b):
        """
        Apply the exponential e^(A dt)*b
        """
        return self._ExpL*b;
    
    def PreCalcSelfFactorization(self):
        """
        Compute the factorization of the operator for inversion. Here it is just 
        a scalar, so we divide by a.
        """
        if (self._Invertible):
            try:
                self._SelfInv = 1/a;   
            except:
                self._Invertible = False;
        
    def addIdentityAndFactor(self,alphaIdentity,alphaL):
        """
        Form the matrix/operator
        A = alphaI*I + alphaL*L
        and precompute its factorization
        """
        self._Id = 1;
        self._SolveOp = alphaIdentity*self._Id+alphaL*self._L;
        self._OpInverse = 1/self._SolveOp;
        
    def SolveSystem(self,b):
        """
        Solve the system
        A x = b, where
        A = alphaI*I + alphaL*L
        has factorization computed above
        """
        return self._OpInverse*b;
    
    def Apply(self,RHS):
        """
        Apply the linear operator L
        """
        return self._L*RHS
        
    def ApplyInverse(self,b):
        """
        Apply L^(-1)
        """
        if (not self._Invertible):
            raise TypeError('Cannot invert a matrix/operator which is not invertible')
        else:
            return self._SelfInv*b;

class DenseLinearOperator(LinearOperator):

    """
    This class implements a dense matrix
    """
    
    def __init__(self,Nunknown):
        """
        Initialize with the number of unknowns. Initialize memory for the matrix representation of L
        """
        self._N = Nunknown;
        self._L = np.zeros((Nunknown,Nunknown));
        
    def setOperator(self,A,invertible=True):
        """
        Set the linear operator A. 
        Here A is a dense matrix
        Invertible = boolean to say whether a is invertible or not
        """
        self._L = A;
        self._Invertible = invertible;
        
    def PreCalcSelfFactorization(self):
        """
        Compute the factorization of the operator for inversion. 
        """
        if (self._Invertible):
            self._LLUFactor = lu_factor(self._L);
        else:
            raise TypeError('Cannot invert a matrix/operator which is not invertible')
    
    def PreCalcExponential(self,dt):
        """
        Precompute the exponential e^(A dt), where dt is 
        the time step size
        """
        self._ExpLDt = expm(self._L*dt);
        
    def ApplyExponential(self,b):
        """
        Apply the exponential e^(A dt)*b
        """
        return np.dot(self._ExpLDt,b);
    
    def addIdentityAndFactor(self,alphaIdentity,alphaL):
        """
        Form the matrix/operator
        A = alphaI*I + alphaL*L
        and precompute its factorization
        
        Accuracy was being lost using Cholesky (only getting 6 digits). Switching
        to LU to get more digits (also makes code simpler)
        """
        self._Id = np.identity(self._N);
        self._MatrixForSolve = alphaIdentity*self._Id+alphaL*self._L;
        self._LUFactor = lu_factor(self._MatrixForSolve);           
    
    def SolveSystem(self,b):
        """
        Solve the system
        A x = b, where
        A = alphaI*I + alphaL*L
        has factorization computed above
        """
        return lu_solve(self._LUFactor,b);
    
    def Apply(self,RHS):
        """
        Apply the linear operator L
        """
        return np.dot(self._L,RHS)
    
    def ApplyInverse(self,b):
        if (not self._Invertible):
            raise TypeError('Cannot invert a matrix/operator which is not invertible')
        else:
            return lu_solve(self._LLUFactor,b); 


class SparseLinearOperator(LinearOperator):

    """
    This class implements a sparse matrix
    """

    def __init__(self,Nunknown):
        """
        Initialize with the number of unknowns. Initialize memory for the matrix representation of L
        """
        self._N = Nunknown;
        self._L = csc_matrix(speye(self._N));
        
    def setOperator(self,A,invertible=True):
        """
        Set the linear operator A. 
        The input here can be a sparse (or dense) matrix in any format. 
        It will be cast to a "csc" matrix,
        which is the required format for LU factorization
        Invertible = boolean to say whether a is invertible or not
        """
        
        self._L = csc_matrix(A); 
        self._Invertible = invertible;
        
    def PreCalcSelfFactorization(self):
        """
        Precompute self factorization (if invertible)
        """
        if (self._Invertible):
            self._LLUFactor = splu(self._L);
        else:
            raise TypeError('Cannot invert a matrix/operator which is not invertible')
        
    def PreCalcExponential(self,dt):
        """
        Precompute the exponential e^(A dt), where dt is 
        the time step size
        """
        self._ExpLDt = sp_expm(self._L*dt);
        
    def ApplyExponential(self,b):
        """
        Apply the exponential e^(A dt)*b
        """
        return self._ExpLDt.dot(b);
        
    def addIdentityAndFactor(self,alphaIdentity,alphaL): 
        """
        Form the matrix/operator
        A = alphaI*I + alphaL*L
        and precompute its factorization
        """
        self._Id = csc_matrix(speye(self._N));
        self._MatrixForSolve = alphaIdentity*self._Id+alphaL*self._L;
        # Cholesky not implemented in scipy for sparse matrices, so resort to LU
        self._LUFactor = splu(self._MatrixForSolve);
        
    def SolveSystem(self,b):
        """
        Solve the system
        A x = b, where
        A = alphaI*I + alphaL*L
        has factorization computed above
        """
        return self._LUFactor.solve(b);
    
    def Apply(self,RHS):
        """
        Apply the linear operator L
        """
        return self._L.dot(RHS);
    
    def ApplyInverse(self,b):
        """
        Apply the inverse of linear operator L
        """
        if (not self._Invertible):
            raise TypeError('Cannot invert a matrix/operator which is not invertible')
        else:
            return self._LLUFactor.solve(b);
            
