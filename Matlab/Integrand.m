classdef Integrand < handle
% Integrand class contains methods for evaluating f(x,t), applying the
% Jacobian of f, solving linear systems with the Jacobian of f, 
% and simple Newton nonlinear systems from implicit ODE integrators

 properties (SetAccess = private)
  evalf         % function handle f(x,t)
  evalJ         % function handle J(x,t)
  n   {mustBeNumeric}    % number of x variables
  Jacobian        % Jacobian df/dx, a LinearOperator object
  tol   {mustBeNumeric}    % maximum relative tolerance required for operations on f or J
  alpha  {mustBeNumeric}    % alpha for solves (default=0), either linear solves (alpha*I+J(x,t))*x=b
     % or nonlinear alpha*x+f(x,t)=b, which uses the linear solve inside Newton's method
  linear_flag {mustBeNumericOrLogical}   % flag for whether f(x,t) is linear in x, i.e. f(x,t) = J(t)*x
 end
 
 methods
 
%------------------------------------------------------------------------------------------
% Main Methods
%------------------------------------------------------------------------------------------

  function obj = Integrand(varargin)
   % Constructor: inputs Evalf, EvalJ, n, J, linear_flag (default = 0), tol (default = 1e-9) 
   Defaults={[],[],[],[],0,1e-9};
   Defaults(1:nargin) = varargin;
   
   % if f(x,t) is NOT linear and no function handle is provided
   if (~strcmpi(class(Defaults{1}),'function_handle')) & (Defaults{5}==0)
    error('Please provide a function handle f(x,t)');
   end
   % if f(x,t) is linear, and no Jacobian is provided
   if (isempty(Defaults{4})) & (Defaults{5}==1)
    error('Please provide a Jacobian J for f(x,t) = J(t)x')
   end
   
   obj.evalf = Defaults{1};
   obj.evalJ = Defaults{2};
   obj.n = Defaults{3};
   obj.Jacobian = Defaults{4}; % in set method, only permits LinearOperator objects 
   obj.linear_flag = Defaults{5};
   obj.tol = Defaults{6};
   obj.alpha = 0.0; % Give default value
   
  end

  function evalf = Evaluate(obj,x,t)
   % Evaluates f=f(x,t) at inputs x, t
   if ~obj.linear_flag
    evalf = obj.evalf(x,t);   % this is assuming f is a function handle
   else
    evalf = obj.ApplyJacobian(x); % this is assuming J is already J=J(x,t). 
    % If J has changed, user must call obj.Jacobian.Update_LinOpA 
   end
  end

  function evalJ = EvaluateJacobian(obj,x,t)
   % Evaluates df/dx=J(x,t) at inputs x, t, and updates the member Jacobian
   if ~obj.linear_flag
    evalJ = obj.evalJ(x,t);   % this is assuming f is a function handle
    obj.Jacobian.Update_LinOpA(evalJ); % tell the LinearOperator J the new value of df/dx
   else
    evalJ = obj.Jacobian.LinOpA; % this is assuming J is already J=J(x,t). 
      % If J has changed, user must call obj.Jacobian.Update_LinOpA 
   end
  end
   
  function Jb = ApplyJacobian(obj,b)
   % Computes J(x,t)*b, taking input b. J is assumed already to be J=J(x,t)
   Jb = obj.Jacobian.ApplyLinOp(b);
  end
  
  function x = SolveLinSys(obj,b,varargin)
   % Solves (alpha*I+J(x,t))*x=b, taking input b and (optional) alpha (otherwise, stored alpha is used)
   if isempty(obj.Jacobian)
    error('Jacobian not provided')
   end
   if nargin > 2
    obj.alpha = varargin{1};
   end
   
   % check if alpha stored in J is the same. If not, update alpha
   if obj.Jacobian.alpha ~= obj.alpha
    obj.Jacobian.Update_alpha(obj.alpha);
   end
   
   x = obj.Jacobian.SolveLinSys(b);
   
  end % SolveLinSys
  
  function x = SolveNonLinSys(obj,b,t,varargin)
   % Solves alpha*x+f(x,t)=b with Newton's method, taking inputs b and t, and 
   % (optional) alpha (otherwise, stored alpha is used),
   % (optional) initial guess x0 (otherwise, x0=0 is used),
   % (optional) max number of Newton iterations max_iter, default=100 
   if nargin > 3
    obj.alpha = varargin{1};
   end

   if obj.Jacobian.alpha ~= obj.alpha
    obj.Jacobian.Update_alpha(obj.alpha);
   end
     
   if obj.linear_flag == 1 % Linear
   
     % The Jacobian object knows how to do this 
     x = obj.SolveLinSys_Jacobian(b);
    
   else % apply Newton's method

     xprev = [];
     if nargin > 4
      xprev = varargin{2};
     end
     if(isempty(xprev))
      xprev = zeros(size(b)); % initial guess is zero
     end 

     kmax = 100;
     if nargin > 5
      kmax = varargin{3};
     end

     % stops either at max iterations, or when absolute tolerance obj.tol is met
     for k=1:kmax
      obj.EvaluateJacobian(x,t);
      xdel = obj.SolveLinSys_Jacobian(-(obj.evalf(xprev,t)+obj.alpha*xprev-b));
      x = xprev+xdel;
      xdelnorm = norm(xdel,inf)
      if xdelnorm < obj.tol
       fprintf('Tolerance reached when k=%d',k);
       break;
      end
      xprev = x;
      if(k==kmax)
       fprintf('Absolute error in Newton solver %g too large after %d iters',xdelnorm, k);
      end
     end % Newton iteration
   
   end % if Linear
   
  end % SolveNonLinSys
  
 end % Methods
 
end % class Integrand

