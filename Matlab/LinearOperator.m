classdef LinearOperator < handle
 % LinearOperator class contains methods for computing b=A*x, solving
 % the linear system (alpha*I+A)*x=b, computing exp(A*t)*x, and an
 % initialize method which precomputes necessary factorizations
 
 properties (Access = public)
  LinOpA  (:,:) {mustBeNumeric}   % n by m matrix operator (can be sparse): R^n -> R^m
  tol  {mustBeNumeric}   % maximum relative tolerance required for operations on A
  alpha  {mustBeNumeric}   % alpha for linear solves (alpha*I+A)*x=b for n=m
  solve_flag {mustBeNumeric}   % user-chosen flag whether to use matinv or LU/Cholesky: 
   % 0 = don't precompute matrices for linear solves, 1 = use Cholesky on alpha*I+A, 
   % -1 = use Cholesky on -(alpha*I+A), other = use LU factorization
   % Note that child classes are free to change the meaning of this flag
  dt   {mustBeNumeric}   % dt for computing exp(A*dt), if dt is set and dt>0
 end
 
 properties (SetAccess = private) 
  expAt  (:,:) {mustBeNumeric}   % n by m matrix exponential of operator
  alphaIA_LU      % decomposition object with the LU factorization of (alpha*I+A) 
  alphaIA_Chol      % R matrix of Cholesky factorization of (alpha*I+A)=RR^T
  alphaIA_Chol_sign     % Is (alpha*I+A) positive definite (+1) or negative definite (-1)
 end
 
 methods
 
%------------------------------------------------------------------------------------------
% Main Methods
%------------------------------------------------------------------------------------------

  function obj = LinearOperator(varargin)
   % Constructor takes inputs: A, solve_flag, tol, alpha, dt
   % alpha and dt can be omitted and set later using the corresponding Update methods
   Defaults={[],0,1e-9,[],[]};
   Defaults(1:nargin) = varargin;
   if isempty(Defaults{1})
    error('Please provide a linear operator')
   end
   obj.LinOpA = Defaults{1};
   obj.solve_flag = Defaults{2};
   obj.tol = Defaults{3};
   obj.alpha = Defaults{4};
   obj.dt = Defaults{5};
   
   % Precompute things to have the object ready for use
   if ~isempty(obj.alpha) % We need to know alpha to initialize
     obj.PrecomputeFactorizations;
   end  
   if ~isempty(obj.dt) % We need to know alpha to initialize
     obj.ComputeExpOp;
   end  

  end
  
  function b = ApplyLinOp(obj,x)
   % Computes b=Ax to tolerance tol. Takes input x.
   % validate x is the correct size
   validateattributes(x,{'double'},{'size',[size(obj.LinOpA,2),1]});
   b = obj.LinOpA*x;
  end

  function b = ApplyExpOp(obj,x)
   % Computes b=exp(A*t)*x for b to tolerance tol. Takes input x.
   % validate x is the correct size
   validateattributes(x,{'double'},{'size',[size(obj.LinOpA,2),1]});
   b = obj.expAt*x;
  end
  
  function x = SolveLinSys(obj,b)
   % Solves (alpha*I+A)*x=b for x to tolerance tol. Takes input b.
   % validate b is the correct size
   validateattributes(b,{'double'},{'size',[size(obj.LinOpA,1),1]});
   
   if obj.solve_flag == 0   % using matlab built-in backslash
    x = (obj.alpha*speye(size(obj.LinOpA,1))+obj.LinOpA)\b;    
   elseif abs(obj.solve_flag) == 1  % using Cholesky factorization
    x = obj.alphaIA_Chol_sign*obj.alphaIA_Chol\(obj.alphaIA_Chol'\b);   
   else      % using LU factorization
    x = obj.alphaIA_LU\b;
   end
  end
  
%------------------------------------------------------------------------------------------
% Update Methods
%------------------------------------------------------------------------------------------ 

  function Update(obj,varargin)
   % Given input matrix A and optional alpha, updates property LinOpA and 
   % precomputes factorizations and matrix exponential if needed.    
   
   if nargin>1
     if(~isempty(varargin{1})); obj.LinOpA = varargin{1}; end
   end
   if nargin>2
    if(~isempty(varargin{2})); obj.alpha = varargin{2}; end
   end   
   if nargin>3
    obj.solve_flag = varargin{2};
   end
   
   obj.PrecomputeFactorizations;
   obj.ComputeExpOp;
  end
  
  function Update_dt(obj,dt)
   % Given input dt>0, updates property dt, precomputing exp(A*dt)
   obj.dt = dt;
   obj.ComputeExpOp;
  end

  function Update_alpha(obj,alpha,varargin)
   % Given input alpha and (optional) solve_flag, updates property alpha for (alpha*I+L)x=b solves 
   % and precomputes factorizations of (alpha*I+L) depending on solve_flag value
   obj.alpha=alpha;
   if nargin>2
    obj.solve_flag = varargin{1};
   end
   obj.PrecomputeFactorizations;
  end  

 end
 
%------------------------------------------------------------------------------------------
% Private, precomputing and computing methods
%------------------------------------------------------------------------------------------

 methods (Access = private)
 
  function PrecomputeFactorizations(obj)
   if obj.solve_flag == 1   % Precompute Cholesky factorization if (alpha*I+A) is SPD
    try
     obj.alphaIA_Chol = chol(obj.alpha*speye(size(obj.LinOpA,1))+obj.LinOpA);
     obj.alphaIA_Chol_sign = +1;
    catch
     error('Cannot compute Cholesky factorization of (alpha*I+A) -- it is not SPD');
    end
   elseif obj.solve_flag == -1   % Precompute Cholesky factorization if -(alpha*I+A) is SPD
    try
     obj.alphaIA_Chol = chol(-(obj.alpha*speye(size(obj.LinOpA,1))+obj.LinOpA));
     obj.alphaIA_Chol_sign = -1;
    catch
     error('Cannot compute Cholesky factorization of -(alpha*I+A) -- it is not SPD');
    end
   elseif obj.solve_flag ~= 0   % Precompute LU factorization as a decompostion object
    obj.alphaIA_LU = decomposition(obj.alpha*speye(size(obj.LinOpA,1))+obj.LinOpA,'lu');
   end
  end
  
  function ComputeExpOp(obj)
   % Computes exp(A*dt) to tolerance tol, if obj.dt exists and is positive 
   if ~isempty(obj.dt)
   if obj.dt > 0
    obj.expAt = expm(obj.LinOpA * obj.dt); % use matlab built-in function for now
   end
   end
  end
  
 end % private methods
  
end % class LinearOperator

